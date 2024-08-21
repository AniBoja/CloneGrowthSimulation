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


# Tree structure for a cell. It points to its children and its parent. The cell id
# is updated at the end of every iteration in the top-level for-loop and reflects
# the order of each evolutionary step
mutable struct cancercell
    mutations::Array{Int64,1}
    fitness::Float64
    epnumber::Float64
    parent::Union{cancercell,Nothing}
    lChild::Union{cancercell,Nothing}
    rChild::Union{cancercell,Nothing}
    status::String
    id::Int64
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
    return cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness),
    copy(cancercellold.epnumber),
    cancercellold.parent, 
    cancercellold.lChild,
    cancercellold.rChild,
    cancercellold.status,
    copy(cancercellold.id))
end

# Start population from a single progenitor cell that has an initial number of mutations, each mutation
# has a unique identifier
function start_population(p, neoep_dist, initial_mut, allowNeoep=true, immThresh=0.5)
    mutID = 1
    N = 1
    cells = cancercell[]
    muts = Dict{Int64,Float64}()
    push!(cells, cancercell([], 1, 0, nothing, nothing, nothing, "extant", 1))
    for i = 1:initial_mut
        if allowNeoep
            cells[1], mutID, neoep_val = newmutations(cells[1], mutID, p, neoep_dist)
        else
            cells[1], mutID, neoep_val = newmutations(cells[1], mutID, 0, neoep_dist)
        end

        muts[mutID-1] = neoep_val
    end

    nonimm = 1 * (cells[1].epnumber < immThresh)

    return cells, mutID, muts, nonimm, N

end

function birthdeath_neoep(b0, d0, Nmax, p, neoep_dist, initial_mut=10, mu=1, immThresh=0.5)

    dmax = d0 #dmax is updated throughout, starts from d0

    #initialize arrays and parameters
    cells, mutID, muts, nonimm, N = start_population(p, neoep_dist, initial_mut)
    Nvec = Int64[]
    push!(Nvec, N)
    t = 0.0
    tvec = Float64[]
    push!(tvec, t)
    nonimmvec = Int64[]
    push!(nonimmvec, nonimm)
    while (N < Nmax) && (t < 300) #set so we can exit simulation where there is a lot of death
        #pick a random cell
        randcell = rand(1:length(cells))
        while cells[randcell].status != "extant"
            randcell = rand(1:length(cells))
        end
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
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - 1 * (cells[randcell].epnumber < immThresh)
            #add new mutations to both new cells, the number of mutations is Poisson distributed
            for i = 1:(rand(Poisson(mu)))
                cells[end-1], mutID, neoep_val = newmutations(cells[end-1], mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            for i = 1:(rand(Poisson(mu)))
                cells[end], mutID, neoep_val = newmutations(cells[end], mutID, p, neoep_dist)
                muts[mutID-1] = neoep_val
            end
            cells[randcell].status = "parent"
            cells[end-1].id = length(cells) - 1
            cells[end-1].parent = cells[randcell]
            cells[end].parent = cells[randcell]
            cells[end].id = length(cells)
            cells[randcell].lChild = cells[end-1]
            cells[randcell].rChild = cells[end]
            
            #note down (non)immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + 1 * (cells[end-1].epnumber < immThresh) + 1 * (cells[end].epnumber < immThresh)

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
            cells[randcell].status = "dead"
            push!(Nvec, N)
            push!(nonimmvec, nonimm)
            Δt = 1 / (Rmax * Nt) .* -log(rand())
            t = t + Δt
            push!(tvec, t)
        end

        #if every cell dies, restart simulation from a single cell again
        if (N == 0)
            cells, mutID, muts, nonimm, N = start_population(p, neoep_dist, initial_mut)
            push!(Nvec, N)
            push!(nonimmvec, nonimm)
            push!(tvec, t)
        end

    end

    return Nvec, tvec, mutID, muts, cells, nonimmvec
end

function process_mutations(cells, detLim)
    mutVec = []
    for i in eachindex(cells)
        append!(mutVec, cells[i].mutations)
    end

    detMutDict = Dict(k => v for (k, v) in countmap(mutVec) if v > detLim)

    println("Mutations processed for ", length(cells), " total cells.")
    return detMutDict

end

function cellnode_to_newick(node::cancercell)
    edge_length = ""
    if (!isnothing(node.parent))
        parent_mutations = Set(node.parent.mutations)
        new_mutations = setdiff(Set(node.mutations), parent_mutations)
        edge_length = length(new_mutations)
    else
        edge_length = 0
    end
    if isnothing(node.lChild) && isnothing(node.rChild)
        return "Cell" * string(node.id) * ":" * string(edge_length)
    else
        lChild_str = !isnothing(node.lChild) ? cellnode_to_newick(node.lChild) : ""
        rChild_str = !isnothing(node.rChild) ? cellnode_to_newick(node.rChild) : ""
        children_str = join(filter(!isempty, [lChild_str, rChild_str]), ",")
        return "($children_str)Cell" * string(node.id) * ":" * string(edge_length)
    end
end

function lineage_to_newick(lineage::Vector{cancercell})
    root_node = lineage[1]
    return cellnode_to_newick(root_node) * ";"
end

# Trim leaves that are either dead or no longer exist
function prune_tree(root::cancercell)
    function prune_recursively!(cell::cancercell)::Bool
        if cell.status == "dead"
            return true
        end
        
        if !isnothing(cell.lChild) && prune_recursively!(cell.lChild)
            cell.lChild = nothing
        end
        
        if !isnothing(cell.rChild) && prune_recursively!(cell.rChild)
            cell.rChild = nothing
        end

        if isnothing(cell.lChild) && isnothing(cell.rChild) && cell.status == "parent"
            return true
        end
        
        return false
    end
    prune_recursively!(root)
end

# Go through each cell -> Identify what mutations it inherited and what mutations it developed -> 
# Write out in a tsv (Cell id, parent mutations, new (or default for progenitor cell) mutations)
function write_tree_mutations(cells, step)
    cell_mutations = open("out/cell_mutations_" * string(step) * ".tsv", "a")
    write(cell_mutations, "id\tparent_mut\tnew_mut\ttip\n")
    write(cell_mutations, string("Cell", cells[1].id, "\t[]", "\t[", join(cells[1].mutations |> collect
                                                                            |> sort, ','), "]", "\tfalse", "\n"))
    for i in eachindex(cells)
        if (!isnothing(cells[i].parent))
            parent_mutations = Set(cells[i].parent.mutations)
            new_mutations = setdiff(Set(cells[i].mutations), parent_mutations)
            tip = cells[i].status == "extant" ? "true" : "false"
            write(cell_mutations, string("Cell", cells[i].id,
                "\t[", join(parent_mutations |> collect |> sort, ','), "]",
                "\t[", join(new_mutations |> collect |> sort, ','), "]",
                "\t",tip,"\n"))
        end
    end
    close(cell_mutations)
end

for i = 1:100
    Nvec, tvec, mutID, muts, cells, immune = birthdeath_neoep(1, d0, popSize, p, neoep_dist, initial_mut, mu)
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    XLSX.writetable("out/preIT_" * string(i) * ".txt", outNDFsim) #Record population size during simulation
    detMutDict = process_mutations(cells, detLim)
    writedlm("out/vaf_preIT_" * string(i) * ".txt", detMutDict) #Save mutation-VAF pairs
    writedlm("out/all_mutations_" * string(i) * ".txt", muts) #Output dictionary storing mutations and their antigenicity
    prune_tree(cells[1])
    write_tree_mutations(cells, i)
    newick_string = lineage_to_newick(cells)
    write("out/newick_" * string(i) * ".tree", newick_string)
end