using Distributions, StatsBase, DataFrames, GLM

type cancercell
    mutations::Array{Int64,1}
    fitness::Int64
end

function newmutations(cancercell, mutID, p)

    numbermutations= 1
    cancercell.mutations = append!(cancercell.mutations, mutID:mutID+numbermutations-1)
    mutID = mutID + numbermutations
    
    neoep = rand()<p
    if neoep
        cancercell.fitness = 0
    end

    return cancercell, mutID, neoep
end

function copycell(cancercellold::cancercell)
  newcancercell::cancercell = cancercell(copy(cancercellold.mutations), copy(cancercellold.fitness))
end

function tumourgrow_birthdeath_neoep(b0, d0, b_, d_, Nmax, p)


    #Rmax starts with b + d and changes once a fitter mutant is introduced, this ensures that
    # b and d have correct units
    Rmax = b0+d_

    #initialize arrays and parameters
    mutID = 1
    cells = cancercell[]
    muts = Int64[]
    neoep_muts = Int64[]
    push!(cells,cancercell([],1))
    for i=1:10
        cells[1],mutID,neoep = newmutations(cells[1],mutID, p)
        if neoep
            push!(neoep_muts, mutID)
        end
        push!(muts,mutID)
    end
    N = 1
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimm = cells[1].fitness
    nonimmvec = Int64[]
    push!(nonimmvec, nonimm)

    while N < Nmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        if cells[randcell].fitness==0 #set death rate according to whether cell is antigenic or not
            d = d_
        else
            d = d0
        end

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < b0

            #population increases by one
            N = N + 1
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - cells[randcell].fitness
            #add new mutations to both new cells
            cells[randcell],mutID,neoep = newmutations(cells[randcell],mutID, p)
            push!(muts,mutID)
            if neoep
                push!(neoep_muts, mutID)
            end
            cells[end],mutID,neoep = newmutations(cells[end],mutID, p)
            push!(muts,mutID)
            if neoep
                push!(neoep_muts, mutID)
            end
            #note down non/immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + cells[randcell].fitness + cells[end].fitness
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        if  (b0+d)<= r
          push!(Nvec, N)
          push!(nonimmvec,nonimm)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - cells[randcell].fitness
            #frequency of cell type decreases
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(Nvec,N)
            push!(nonimmvec, nonimm)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies reinitialize simulation
        if (N == 0)
                mutID = 1
        cells = cancercell[]
        muts = Int64[]
    neoep_muts = Int64[]
    push!(cells,cancercell([],1))
            for i=1:10
                cells[1],mutID,neoep = newmutations(cells[1],mutID, p)
                if neoep
                    push!(neoep_muts, mutID)
                end
            end
    N = 1
    push!(Nvec,N)
        nonimm = cells[1].fitness
    push!(nonimmvec, nonimm)
    push!(tvec,t)
        end

    end
    
    return Nvec, tvec, mutID, neoep_muts, cells, nonimmvec
end

function tumourgrow_post_it(b0, d0, b_, d_, Tmax, startPopSize, cells, mutID, neoep_muts, fitnessChange)
    #Function only to track the size of population afterwards, don't care about mutants
    #Either immunogenicity status is fixed and heritable or new mutations can change it
    #Follow through which mutations are neo-epitopes from previous simulation

    #Rmax starts with b + d
    Rmax = b0+d_

    #initialize arrays and parameters: cells are inputed, only time and pop-size matters
    startPopSize = length(cells)
    N = startPopSize
    Nvec = Int64[]
    push!(Nvec,N)
    t = 0.0
    tvec = Float64[]
    push!(tvec,t)
    nonimm = 0
    for i=1:length(cells)
        nonimm = nonimm + cells[i].fitness
    end
    nonimmvec = Int64[]
    push!(nonimmvec, nonimm)

    while t < Tmax

        #pick a random cell
        randcell = rand(1:N)
        r = rand(Uniform(0,Rmax))
        Nt = N
        
        if cells[randcell].fitness==0 #set death rate according to whether cell is antigenic or not
            d = d_
        else
            d = d0
        end

        #birth event if r<birthrate, access correct birthrate from cells array
        if r < b0

            #population increases by one
            N = N + 1
            #total fitness (nonimmunogenicity) decreases as it might change in mutation step
            nonimm = nonimm - cells[randcell].fitness
            #copy cell and mutations for cell that reproduces
            push!(cells, copycell(cells[randcell]))
            if fitnessChange == 1
                #add new mutations to both new cells but don't note mutations anymore (nevertheless cell can change fitness state)
                cells[randcell],mutID,neoep = newmutations(cells[randcell],mutID, p)
                if neoep
                    push!(neoep_muts, mutID)
                end
                cells[end],mutID,neoep = newmutations(cells[end],mutID, p)
                if neoep
                    push!(neoep_muts, mutID)
                end
            end
            #note down non/immunogenicity stored in fitness for the new cells:
            nonimm = nonimm + cells[randcell].fitness + cells[end].fitness
            
            push!(nonimmvec, nonimm)
            push!(Nvec, N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
            
        end

        if  (b0+d)<= r
          push!(nonimmvec, nonimm)
          push!(Nvec, N)
          Δt =  1/(Rmax * Nt) .* - log(rand())
          t = t + Δt
          push!(tvec,t)
        end

        #death event if b<r<b+d
        if b0 <= r < (b0+d)

            #population decreases by 1, overall fitness score also decreases if it was non-zero
            N = N - 1
            nonimm = nonimm - cells[randcell].fitness
            #remove deleted cell
            deleteat!(cells,randcell)
            push!(nonimmvec, nonimm)
            push!(Nvec,N)
            Δt =  1/(Rmax * Nt) .* - log(rand())
            t = t + Δt
            push!(tvec,t)
        end

        #every cell dies, automatically jump to Tmax
        if (N == 0)
          t = Tmax
        end

       #if we reach a high population showing exponential growth, also jump to Tmax
       if (N > 1.5*startPopSize)
         t = Tmax
       end

    end
    
    return Nvec, tvec, cells,nonimmvec, neoep_muts, mutID
end

function process_mutations(cells, neoep_muts, detLim, popSize)
	fitnessVec = []
	mutVec = []
    for i=1:length(cells)
        append!(mutVec, cells[i].mutations)
    end

    detMutDict = filter((k, v) -> v > detLim, countmap(mutVec))
    detEpMutDict = filter((k, v) -> k in neoep_muts, detMutDict)
    VAF = collect(values(detMutDict))/(2*popSize)
    VAFtotal = Array{Float64}(0)
    append!(VAFtotal, VAF)

    VAFep = collect(values(detEpMutDict))/(2*popSize)
    println("Mutations processed for ", length(cells), "cells.")
	return detMutDict, detEpMutDict, VAFtotal, VAFep

end

function analyse_vaf(VAFtotal, steps, fmin, fmax)
	steps = fmax:-0.0001:fmin
	cumcount = Array{Int64}(0)
	cumfm = Array{Float64}(0)
	invv = Array{Float64}(0)
	for i in steps
	    push!(cumcount, sum(VAFtotal .>= i))
	    push!(cumfm, sum(VAFtotal[VAFtotal .>= i]))
	    push!(invv, 1/i-1/fmax)
	end
	cumcount = cumcount-cumcount[1]
	cumfm = cumfm - cumfm[1]
	vafDF = DataFrame(invf = invv, cumcount = cumcount)
	lmfit = lm(@formula(cumcount ~ invf+0), vafDF)
	rsq = r2(lmfit)

	return cumcount, cumfm, rsq
end



outRvec = []
outREpvec = []
outMutRvec = []
outvafDF = DataFrame(invf = (1./steps - 1/fmax))
outvafEpDF = DataFrame(invf = (1./steps - 1/fmax))
outfimDF = DataFrame(f=steps)
outfimEpDF = DataFrame(f=steps)
for i=1:200
    Nvec, tvec, mutID, neoep_muts, cells, immune = tumourgrow_birthdeath_neoep(1, d0, 1, d_, popSize, p);
    outNDFsim = DataFrame(t=tvec, N=Nvec, nonImm=immune)
    writetable("preIT_"*string(i)*".txt", outNDFsim) #Record population size during simulation

    detMutDict, detEpMutDict, VAFtotal, VAFep = process_mutations(cells, neoep_muts, detLim, popSize)
    mutR = length(detEpMutDict)/length(detMutDict)
    writedlm("vaf_preIT_"*string(i)*".txt",detMutDict) #Save mutation-VAF pairs
    #Process mutations, i.e. cumulative distributions
    cumcount, cumfm, rsq = analyse_vaf(VAFtotal, steps, fmin, fmax)
    cumcountEp, cumfmEp, rsqEp = analyse_vaf(VAFep, steps, fmin, fmax)
    push!(outMutRvec, mutR)
    push!(outRvec, rsq)
    push!(outREpvec, rsqEp)
    outvafDF[Symbol("sim_$(i)")] = cumcount
    outfimDF[Symbol("sim_$(i)")] = cumfm
    outvafEpDF[Symbol("sim_$(i)")] = cumcountEp
    outfimEpDF[Symbol("sim_$(i)")] = cumfmEp

    #After all is recorded, we do 'generic' therapy
    NvecT, tvecT, cellsT, immuneT, neoep_mutsT, mutIDT = tumourgrow_post_it(1, d_t, 1, d_t, TmaxT, popSize, cells, mutID, neoep_muts, 1)
    outTDFsim = DataFrame(t=tvecT, N=NvecT, nonImm=immuneT )
    writetable("midIT_"*string(i)*".txt", outTDFsim)
    #Process the mutation after therapy just like before
    popSizeT = length(cellsT) 
    detMutDictT, detEpMutDictT, VAFtotalT, VAFepT = process_mutations(cellsT, neoep_mutsT, detLim, popSizeT)
    writedlm("vaf_midIT_"*string(i)*".txt",detMutDictT)

    #Then we go on to simulate the population after 'immunotherapy' and process mutations
    NvecIT, tvecIT, cellsIT, immuneIT, neoep_mutsIT = tumourgrow_post_it(1, d0_it, 1, d_it, TmaxIT, popSizeT, cellsT, mutIDT, neoep_mutsT, 1)
    outITDFsim = DataFrame(t=tvecIT, N=NvecIT, nonImm=immuneIT )
    writetable("postIT_"*string(i)*".txt", outITDFsim)
    popSizeIT = length(cellsIT)
    detMutDictIT, detEpMutDictIT, VAFtotalIT, VAFepIT = process_mutations(cellsIT, neoep_mutsIT, detLim, popSizeIT)
    writedlm("vaf_postIT_"*string(i)*".txt",detMutDictIT) #Save mutation-VAF pairs post IT
    writedlm("neoep_mutations_"*string(i)*".txt", neoep_mutsIT)

end

writetable("Vafdf.csv",outvafDF)
writetable("Fimdf.csv",outfimDF)
writetable("Vafdf_ep.csv",outvafEpDF)
writetable("Fimdf_ep.csv",outfimEpDF)
writedlm("Rsqs.txt",outRvec)
writedlm("Rsqs_ep.txt", outREpvec)
writedlm("MutRatios.txt", outMutRvec)
