using DADA2, Random

const SCIENCE_ROOT = joinpath(@__DIR__,"generated","science")
mkpath(SCIENCE_ROOT)
const SCENARIOS = ["no_errors","calibrated","rare_snv","miscalibrated","indels","chimeras"]
const SEEDS = [11,23,37,51,71]

function simulator_errors(scale)
    # Independently constructed asymmetric transition model, not phredErrors().
    rates = zeros(16,41)
    for q in 0:40, source in 1:4
        error = min(0.25,scale*10.0^(-q/10))
        alternatives = [b for b in 1:4 if b != source]
        weights = Float64[1,2,3]
        rates[4(source-1)+source,q+1] = 1-error
        for (j,b) in enumerate(alternatives)
            rates[4(source-1)+b,q+1] = error*weights[j]/sum(weights)
        end
    end
    rates
end
function mutate_read(rng,sequence,quality,rates; indels=false)
    out = UInt8[]; qual = Int[]
    bases = UInt8.(['A','C','G','T'])
    for (pos,nt) in enumerate(codeunits(sequence))
        # Independent Bernoulli insertion/deletion process for stress tests.
        indels && rand(rng) < 0.001 && continue
        source = findfirst(==(nt),bases)
        u = rand(rng); cumulative = 0.0; observed = 4
        for b in 1:4
            cumulative += rates[4(source-1)+b,quality[pos]+1]
            if u <= cumulative
                observed = b; break
            end
        end
        push!(out,bases[observed]); push!(qual,quality[pos])
        if indels && rand(rng) < 0.001
            push!(out,rand(rng,bases)); push!(qual,quality[pos])
        end
    end
    String(out),qual
end
function changed_at(sequence,pos)
    bytes = collect(codeunits(sequence))
    bytes[pos] = bytes[pos] == 0x41 ? 0x43 : 0x41
    String(bytes)
end

open(joinpath(SCIENCE_ROOT,"manifest.tsv"),"w") do manifest
    println(manifest,"case\tscenario\tseed\tbiological_reads\tartifact_reads")
    for scenario in SCENARIOS, seed in SEEDS
        rng = MersenneTwister(seed)
        base = String(rand(rng,UInt8.(['A','C','G','T']),150))
        truth = [base,String(rand(rng,UInt8.(['A','C','G','T']),150)),
                 String(rand(rng,UInt8.(['A','C','G','T']),150)),changed_at(base,60),changed_at(base,100)]
        counts = scenario == "rare_snv" ? [5000,1000,200,10,1] : [3000,1000,500,100,20]
        if scenario == "indels"
            truth[5] = base[1:90]*"A"*base[91:end] # genuine length variant, plus noisy indels
        end
        scale = scenario == "no_errors" ? 0.0 : scenario == "miscalibrated" ? 5.0 : 1.0
        rates = simulator_errors(scale)
        name = "$(scenario)_$(seed)"; dir = joinpath(SCIENCE_ROOT,name); mkpath(dir)
        open(joinpath(dir,"truth.tsv"),"w") do io
            println(io,"sequence\tcount")
            for (s,n) in zip(truth,counts)
                println(io,s,'\t',n)
            end
        end
        open(joinpath(dir,"oracle.tsv"),"w") do io
            for row in eachrow(rates)
                println(io,join(row,'\t'))
            end
        end
        templates = copy(truth); totals = copy(counts)
        artifact_reads = scenario == "chimeras" ? 200 : 0
        if artifact_reads > 0
            chimera = truth[1][1:75]*truth[2][76:end]
            push!(templates,chimera); push!(totals,artifact_reads)
            write(joinpath(dir,"artifacts.tsv"),"sequence\tcount\n$chimera\t$artifact_reads\n")
        end
        open(joinpath(dir,"reads.fastq"),"w") do io
            readid = 0
            for (sequence,n) in zip(templates,totals), _ in 1:n
                quality = rand(rng,25:40,length(sequence))
                read,qs = mutate_read(rng,sequence,quality,rates;indels=scenario=="indels")
                readid += 1
                println(io,"@read",readid,"\n",read,"\n+")
                println(io,String(UInt8.(qs .+ 33)))
            end
        end
        println(manifest,name,'\t',scenario,'\t',seed,'\t',sum(counts),'\t',artifact_reads)
    end
end
println("Generated ",length(SCENARIOS)*length(SEEDS)," simulated conditions across five seeded communities at ",SCIENCE_ROOT)
