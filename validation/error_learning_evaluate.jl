using DADA2, Statistics

root = joinpath(@__DIR__,"generated","science")
modelroot = joinpath(@__DIR__,"generated","error_learning")
mkpath(modelroot)
function readcounts(path)
    Dict(s=>parse(Int,n) for (s,n) in split.(readlines(path)[2:end],'\t'))
end
function matrixfile(path)
    permutedims(hcat([parse.(Float64,split(line,'\t')) for line in readlines(path)]...))
end
function baseweights(truth)
    bases = "ACGT"
    totals = [sum(n*count(==(base),s) for (s,n) in truth) for base in bases]
    totals ./ sum(totals)
end
function calibration(rates,oracle,weights)
    # The generator assigns Q25:40 uniformly, independently of the template.
    # Do not weight by inferred transition exposures: these may be biased.
    absolute = 0.0; predicted = 0.0; actual = 0.0; squaredlog = Float64[]
    for from in 1:4, to in 1:4
        from == to && continue
        row = 4(from-1)+to
        for q in 25:40
            p = rates[row,q+1]; truth = oracle[row,q+1]; w = weights[from]/16
            absolute += w*abs(p-truth)
            predicted += w*p; actual += w*truth
            truth > 0 && push!(squaredlog,(log10(p)-log10(truth))^2)
        end
    end
    (absolute, isempty(squaredlog) ? "NA" : sqrt(mean(squaredlog)),
     predicted, actual, actual == 0 ? "NA" : predicted/actual)
end
function metrics(truth,actual)
    tp = length(intersect(keys(truth),keys(actual)))
    fp = length(setdiff(keys(actual),keys(truth)))
    eligible = [s for s in keys(truth) if truth[s] >= 2]
    precision = isempty(actual) ? 0.0 : tp/length(actual)
    recall = count(s -> haskey(actual,s),eligible)/length(eligible)
    ntruth = sum(values(truth)); nactual = sum(values(actual))
    tv = nactual == 0 ? 1.0 : 0.5sum(abs(get(truth,s,0)/ntruth-get(actual,s,0)/nactual) for s in union(keys(truth),keys(actual)))
    (tp,fp,precision,recall,tv,nactual)
end
cases = [split(line,'\t') for line in readlines(joinpath(root,"manifest.tsv"))[2:end]]
data = Dict(row[1]=>derepFastq(joinpath(root,row[1],"reads.fastq");keepmap=false) for row in cases)
output = joinpath(@__DIR__,"error_learning_results.tsv")
open(output,"w") do io
    println(io,"case\tscenario\tseed\tmethod\tsubstitution_l1\tlog10_rmse\tpredicted_error_rate\ttrue_error_rate\terror_ratio\ttp\tfp\tprecision\trecall_ge2\tabundance_tv\taccepted_reads\terror_converged\tinference_converged\titerations")
    function evaluate(row,method,err;learned=nothing,reference=false)
        name,scenario,seed = row[1:3]
        dir = joinpath(root,name); truth = readcounts(joinpath(dir,"truth.tsv"))
        oracle = matrixfile(joinpath(dir,"oracle.tsv"))
        if learned !== nothing
            open(joinpath(modelroot,"$(name)_$(method).tsv"),"w") do model
                for values in eachrow(err.rates)
                    println(model,join(values,'\t'))
                end
            end
        end
        if reference
            actual = readcounts(joinpath(dir,"r_learned.tsv")); inference = "NA"
        else
            result = dada(data[name],err)
            table = removeBimeraDenovo(makeSequenceTable([result]);method=:pooled)
            actual = Dict(s=>n for (s,n) in zip(table.sequences,vec(table.counts)))
            inference = result.converged
        end
        values = (name,scenario,seed,method,calibration(err.rates,oracle,baseweights(truth))...,
            metrics(truth,actual)...,learned === nothing ? "NA" : learned.converged,
            inference, learned === nothing ? "NA" : learned.iterations)
        println(io,join(values,'\t')); flush(io)
    end
    for row in cases
        println("Error calibration: ",row[1]); flush(stdout)
        for smoother in (:linear,:loess)
            learned = learnErrors(data[row[1]];smoother,maxIter=10)
            evaluate(row,"julia_$(smoother)",learned.errors;learned)
        end
        learned = learnErrors(joinpath(root,row[1],"reads.fastq");maxIter=10)
        evaluate(row,"julia_loess_raw",learned.errors;learned)
        evaluate(row,"r_learned",ErrorModel(matrixfile(joinpath(root,row[1],"r_errors.tsv")));reference=true)
    end
    for scenario in unique(row[2] for row in cases)
        group = filter(row -> row[2] == scenario,cases)
        println("Shared error calibration: ",scenario); flush(stdout)
        for smoother in (:linear,:loess)
            learned = learnErrors([data[row[1]] for row in group];smoother,maxIter=10)
            for row in group
                evaluate(row,"julia_$(smoother)_shared",learned.errors;learned)
            end
        end
        learned = learnErrors([joinpath(root,row[1],"reads.fastq") for row in group];maxIter=10)
        for row in group
            evaluate(row,"julia_loess_raw_shared",learned.errors;learned)
        end
    end
end
println("Saved ",output)
