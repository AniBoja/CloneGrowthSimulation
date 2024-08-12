#=
    Julia code implementing simulations of a stochastic branching process model of tumour growth and mutation accumulations, under negative selection from the immune system due to randomly arising antigenic mutations. This code corresponds to an ideal population/measurement, with no measurement noise (but intrinsic stochasticity) and user-defined detection limit for mutations.

    For details of the model see Lakatos et al, biorXiv, 2019: https://www.biorxiv.org/content/10.1101/536433v1 and the Readme at https://github.com/elakatos/CloneGrowthSimulation

    This file implements the base simulations without immune escape and variable mutation rate and selection strength.
    Parameters: d0, popSize, detLim, p, initial_mut, mu, neoep_dist
    Outputs:
        - all_mutations_<i>.txt : every mutation and its antigenicity value in two tab-separated columns
        - preIT_<i>.txt : cell growth curve with time, total population and non-immunogenic population in three comma-separated columns
        - vaf_preIT_<i>.txt : every mutation present in >detLim cells and the number of cells its present, two tab-separated columns

    @Author: Eszter Lakatos (e.lakatos@qmul.ac.uk)
    Inspired by Marc J. Williams' CancerSeqSim (https://github.com/marcjwilliams1/CancerSeqSim.jl)
=#

using Distributions, StatsBase, DataFrames, GLM, XLSX, DelimitedFiles

function getFitness(n)
    (1 + s * n)
end



mutable struct cancercell
    mutations::Array{Int64,1}
    fitness::Float64
    epnumber::Float64
end

# Tree structure for a cell. It points to its children and its parent. The cell id
# is updated at the end of every iteration in the top-level for-loop and reflects
# the order of each evolutionary step
mutable struct cellnode
    parent::Union{cellnode,Nothing}
    lChild::Union{cellnode,Nothing}
    rChild::Union{cellnode,Nothing}
    cellid::Union{Int64,String}
    cell::cancercell
end

mutable struct lineage
    cancercells::Vector{cellnode}
end

function newmutations(cancercell, mutID, p, neoep_dist)
    cancercell.mutations = append!(cancercell.mutations, mutID)
    mutID = mutID + 1

    neoep = rand() < p
    if neoep
        neoep_value = max.(0, rand(neoep_dist)) #Take distribution from top-defined params
        cancercell.epnumber = cancercell.epnumber + neoep_value
        cancercell.fitness = getFitness(cancercell.epnumber) #fitness is affected by the number of mutations
    else
        neoep_value = 0
    end

    return cancercell, mutID, neoep_value
end

function copycell(cancercellold::cancercell)
    newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness), copy(cancercellold.epnumber))
end

# Start population from a single progenitor cell that has an initial number of mutations, each mutation
# has a unique identifier
function start_population(p, neoep_dist, initial_mut, allowNeoep=true, immThresh=0.5)
    mutID = 1
    N = 1
    cells = cancercell[]
    muts = Dict{Int64,Float64}()
    push!(cells, cancercell([], 1, 0))
    cell = cellnode(nothing, nothing, nothing, 1, cancercell([], 1, 0))
    history = lineage([cell])
    for i = 1:initial_mut
        if allowNeoep
            cells[1], mutID, neoep_val = newmutations(cells[1], mutID, p, neoep_dist)
        else
            cells[1], mutID, neoep_val = newmutations(cells[1], mutID, 0, neoep_dist)
        end

        muts[mutID-1] = neoep_val
    end

    nonimm = 1 * (cells[1].epnumber < immThresh)

    return cells, mutID, muts, nonimm, N, history

end

function birthdeath_neoep(b0, d0, Nmax, p, neoep_dist, initial_mut=10, mu=1, immThresh=0.5)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, muts, nonimm, N, history = start_population(p, neoep_dist, initial_mut)
    Nvec = Int64[]
    push!(Nvec, N)
    t = 0.0
    tvec = Float64[]
    push!(tvec, t)
    nonimmvec = Int64[]
    push!(nonimmvec, nonimm)
    while (N < Nmax) & (t < 300) #set so we can exit simulation where there is a lot of death
        #pick a random cell
        randcell = rand(1:N)
        Nt = N

        #a cell's immunogenicity depends on its fitness, i.e. the summed antigenicity of neoepitopes
        d = max(0, (d0 - b0) * cells[randcell].fitness + b0)

        if (d > dmax) #update dmax to keep track of the highest death rate in the whole population
            dmax = d
        end

        Rmax = b0 + dmax

        r = rand(Uniform(0, Rmax)) #Pick which reaction should happen to cell     

        # If r < birthrate, a birth event happens: a new cell is created and randcell updated as a new one
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1 * (cells[randcell].epnumber < immThresh)

            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i = 1:(rand(Poisson(mu)))
                cells[randcell], mutID, neoep_val = newmutations(cells[randcell], mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            for i = 1:(rand(Poisson(mu)))
                cells[end], mutID, neoep_val = newmutations(cells[end], mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            for i in eachindex(history.cancercells)
                # The cell ids are temporarily assigned the cell's position in the array "cells"
                if (history.cancercells[i].cellid == randcell)
                    lChild = cellnode(history.cancercells[i],
                        nothing,
                        nothing,
                        randcell,
                        cells[randcell])
                    rChild = cellnode(history.cancercells[i],
                        nothing,
                        nothing,
                        length(cells),
                        cells[end])
                    history.cancercells[i].lChild = lChild
                    history.cancercells[i].rChild = rChild
                    push!(history.cancercells, lChild, rChild)
                    # Once a cell has divided, the original cell can be ignored and is marked as "parent"
                    history.cancercells[i].cellid = "parent"
                    break
                end
            end

            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1 * (cells[randcell].epnumber < immThresh) + 1 * (cells[end].epnumber < immThresh)

            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt = 1 / (Rmax * Nt) .* -log(rand())
            t = t + Δt
            push!(tvec, t)

        end

        #if r has neither birth or death (only possible if it is a non-immunogenic cell), nothing happens
        if (b0 + d) <= r
            push!(Nvec, N)
            push!(nonimmvec, nonimm)
            Δt = 1 / (Rmax * Nt) .* -log(rand())
            t = t + Δt
            push!(tvec, t)
        end

        #death event if r > b but < d
        if b0 <= r < (b0 + d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - 1 * (cells[randcell].epnumber < immThresh)
            #remove deleted cell
            deleteat!(cells, randcell)
            push!(Nvec, N)
            push!(nonimmvec, nonimm)
            Δt = 1 / (Rmax * Nt) .* -log(rand())
            t = t + Δt
            push!(tvec, t)
        end

        #if every cell dies, restart simulation from a single cell again
        if (N == 0)
            cells, mutID, muts, nonimm, N, history = start_population(p, neoep_dist, initial_mut)
            push!(Nvec, N)
            push!(nonimmvec, nonimm)
            push!(tvec, t)
        end

    end

    return Nvec, tvec, mutID, muts, cells, nonimmvec, history
end

function process_mutations(cells, detLim)
    mutVec = []
    for i in eachindex(cells)
        append!(mutVec, cells[i].mutations)
    end

    detMutDict = Dict(k => v for (k, v) in countmap(mutVec) if v > detLim)

    println("Mutations processed for ", length(cells), " cells.")
    return detMutDict

end

function cellnode_to_newick(node::cellnode)
    if isnothing(node.lChild) && isnothing(node.rChild)
        return string(node.cellid)
    else
        lChild_str = !isnothing(node.lChild) ? cellnode_to_newick(node.lChild) : ""
        rChild_str = !isnothing(node.rChild) ? cellnode_to_newick(node.rChild) : ""
        children_str = join(filter(!isempty, [lChild_str, rChild_str]), ",")
        return "($children_str)" * string(node.cellid)
    end
end

function lineage_to_newick(history::lineage)
    root_node = history.cancercells[1]
    return cellnode_to_newick(root_node) * ";"
end

# Go through each cell -> Identify what mutations it inherited and what mutations it developed -> 
# Write out in a csv (Cell id, parent mutations, new (or default for progenitor cell) mutations)
function write_tree_mutations(history, step)
    cell_count = 0
    cell_mutations = open("out/cell_mutations" * string(step) * ".csv", "a")
    write(cell_mutations, "id,parent_mut,child_mut\n")
    write(cell_mutations, string(0, ",[]", ",[", join(history.cancercells[1].cell.mutations |> collect
                                                                            |> sort, ' '), "]", "\n"))
    for j in eachindex(history.cancercells)
        history.cancercells[j].cellid = cell_count
        cell_count += 1
        if (!isnothing(history.cancercells[j].parent))
            parent_mutations = Set(history.cancercells[j].parent.cell.mutations)
            new_mutations = setdiff(Set(history.cancercells[j].cell.mutations), parent_mutations)
            write(cell_mutations, string(history.cancercells[j].cellid,
                ",[", join(parent_mutations |> collect |> sort, ' '), "]",
                ",[", join(new_mutations |> collect |> sort, ' '), "]\n"))
        end
    end
    close(cell_mutations)
end

for i = 1:100
    Nvec, tvec, mutID, muts, cells, immune, history = birthdeath_neoep(1, d0, popSize, p, neoep_dist, initial_mut, mu)
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    XLSX.writetable("out/preIT_" * string(i) * ".txt", outNDFsim) #Record population size during simulation
    detMutDict = process_mutations(cells, detLim)
    writedlm("out/vaf_preIT_" * string(i) * ".txt", detMutDict) #Save mutation-VAF pairs
    writedlm("out/all_mutations_" * string(i) * ".txt", muts) #Output dictionary storing mutations and their antigenicity
    write_tree_mutations(history, i)
    newick_string = lineage_to_newick(history)
    write("out/newick_" * string(i) * ".tree", newick_string)
end