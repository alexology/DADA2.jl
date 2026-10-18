using DADA2

root = joinpath(@__DIR__,"generated","science")
output = joinpath(@__DIR__,"scientific_results.tsv")
function read_counts(path)
    result = Dict{String,Int}()
    for line in readlines(path)[2:end]
        s,n = split(line,'\t'); result[s] = parse(Int,n)
    end
    result
end
function read_matrix(path)
    rows = [parse.(Float64,split(line,'\t')) for line in readlines(path)]
    ErrorModel(permutedims(hcat(rows...)))
end
function write_counts(path,table)
    open(path,"w") do io
        println(io,"sequence\tcount")
        for j in eachindex(table.sequences)
            println(io,table.sequences[j],'\t',table.counts[1,j])
        end
    end
end
function metrics(truth,actual)
    tp = length(intersect(keys(truth),keys(actual)))
    fp = length(setdiff(keys(actual),keys(truth)))
    fn = length(truth)-tp
    precision = tp+fp == 0 ? 0.0 : tp/(tp+fp)
    recall = tp/length(truth)
    eligible = [s for s in keys(truth) if truth[s] >= 2]
    recall_ge2 = count(s -> haskey(actual,s),eligible)/length(eligible)
    ntruth = sum(values(truth)); nactual = sum(values(actual))
    universe = union(keys(truth),keys(actual))
    tv = nactual == 0 ? 1.0 : 0.5sum(abs(get(truth,s,0)/ntruth-get(actual,s,0)/nactual) for s in universe)
    spurious = nactual == 0 ? 0.0 : sum(get(actual,s,0) for s in setdiff(keys(actual),keys(truth));init=0)/nactual
    (tp,fp,fn,precision,recall,recall_ge2,tv,spurious,nactual)
end
rows = readlines(joinpath(root,"manifest.tsv"))[2:end]
open(output,"w") do out
    println(out,"case\tscenario\tseed\tmethod\ttp\tfp\tfn\tprecision\trecall\trecall_ge2\tabundance_tv\tspurious_read_fraction\taccepted_reads\tinference_converged\terror_converged")
    for row in rows
        name,scenario,seed,_,_ = split(row,'\t')
        println("Julia validation: ",name)
        dir = joinpath(root,name); truth = read_counts(joinpath(dir,"truth.tsv"))
        d = derepFastq(joinpath(dir,"reads.fastq");keepmap=false)
        for mode in ("oracle","learned","r_errors")
            error_converged = "NA"
            if mode == "learned"
                learned = learnErrors(d;maxIter=10,threaded=true)
                err = learned.errors; error_converged = string(learned.converged)
            else
                path = joinpath(dir,mode=="oracle" ? "oracle.tsv" : "r_errors.tsv")
                isfile(path) || error("Missing reference error model: $path; inspect r_status.tsv")
                err = read_matrix(path)
            end
            result = dada(d,err;threaded=true)
            table = removeBimeraDenovo(makeSequenceTable([result]);method=:pooled)
            path = joinpath(dir,"julia_$mode.tsv"); write_counts(path,table)
            actual = read_counts(path)
            println(out,join((name,scenario,seed,"julia_$mode",metrics(truth,actual)...,result.converged,error_converged),'\t'))
        end
        for mode in ("oracle","learned")
            path = joinpath(dir,"r_$mode.tsv"); isfile(path) || error("Missing R result: $path; inspect r_status.tsv")
            actual = read_counts(path)
            println(out,join((name,scenario,seed,"r_$mode",metrics(truth,actual)...,"NA","NA"),'\t'))
        end
    end
end
println("Scientific metrics: ",output)
