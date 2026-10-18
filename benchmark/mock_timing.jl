const SCRIPT_START = time_ns()
using DADA2, LinearAlgebra
BLAS.set_num_threads(1)
const LOAD_SECONDS = (time_ns()-SCRIPT_START)/1e9

function readmodel(path)
    ErrorModel(permutedims(hcat([parse.(Float64,split(line,'\t')) for line in readlines(path)]...)))
end
function timedstep(f)
    started = time_ns()
    result = f()
    result,(time_ns()-started)/1e9
end
function pipeline(root,dir,mode,shared)
    f = joinpath(dir,"F.fastq.gz"); r = joinpath(dir,"R.fastq.gz")
    stats,filtering = timedstep(() -> filterAndTrim(
        joinpath(root,"Mock_S280_L001_R1_001.fastq"),joinpath(root,"Mock_S280_L001_R2_001.fastq"),f,r;
        optionsF=FilterOptions(truncLen=240,maxEE=2.0),optionsR=FilterOptions(truncLen=160,maxEE=2.0)))
    dereps,dereplication = timedstep(() -> (derepFastq(f),derepFastq(r)))
    fd,rd = dereps
    iterationsF = iterationsR = 0; converged = true
    if mode == "shared_errors"
        ef,er = shared; learning = 0.0
    else
        learned,learning = timedstep(() -> mode == "full_raw" ?
            (learnErrors(f;maxIter=10,omegaC=0),learnErrors(r;maxIter=10,omegaC=0)) :
            (learnErrors(fd;maxIter=10,omegaC=0),learnErrors(rd;maxIter=10,omegaC=0)))
        lf,lr = learned; ef = lf.errors; er = lr.errors
        iterationsF = lf.iterations; iterationsR = lr.iterations
        converged = lf.converged && lr.converged
    end
    inferred,inference = timedstep(() -> (dada(fd,ef),dada(rd,er)))
    fr,rr = inferred
    merged,merging = timedstep(() -> mergePairs(fr,rr;minOverlap=12,maxMismatch=0))
    clean,community = timedstep() do
        table = makeSequenceTable([merged])
        keep = [250 <= length(s) <= 256 for s in table.sequences]
        table = SequenceTable(table.samples,table.sequences[keep],table.counts[:,keep])
        removeBimeraDenovo(table;method=:pooled,minFoldParentOverAbundance=1.0)
    end
    (stats=stats,table=clean,times=(filtering,dereplication,learning,inference,merging,community),
     iterations=(iterationsF,iterationsR),converged=converged && fr.converged && rr.converged)
end
function main()
    length(ARGS) == 1 || error("Pass the output TSV path")
    root = joinpath(@__DIR__,"..","validation","generated","public")
    countsdir = joinpath(@__DIR__,"..","validation","generated","mock_timing")
    mkpath(countsdir)
    shared = (readmodel(joinpath(root,"mock_errF.tsv")),readmodel(joinpath(root,"mock_errR.tsv")))
    open(ARGS[1],"w") do io
        println(io,"implementation\tthreads\tmode\trepetition\tphase\tfilter_seconds\tderep_seconds\tlearning_seconds\tinference_seconds\tmerge_seconds\tcommunity_seconds\ttotal_seconds\tscript_load_seconds\tfirst_script_seconds\treads_in\tfiltered\tasvs\tfinal_reads\terror_iterations_F\terror_iterations_R\tconverged")
        for mode in ("full_mean","shared_errors","full_raw"), rep in 0:5
            GC.gc() # Outside the measured pipeline, as with gc() in R.
            result,total = mktempdir() do dir
                timedstep(() -> pipeline(root,dir,mode,shared))
            end
            firstseconds = mode == "full_mean" && rep == 0 ? (time_ns()-SCRIPT_START)/1e9 : "NA"
            table = result.table
            @assert (result.stats.reads_in,result.stats.reads_out,length(table.sequences),sum(table.counts)) == (4779,4314,20,4269)
            @assert result.converged
            values = ("Julia",Threads.nthreads(),mode,rep,rep == 0 ? "first_call" : "warm",
                result.times...,total,LOAD_SECONDS,firstseconds,result.stats.reads_in,result.stats.reads_out,
                length(table.sequences),sum(table.counts),result.iterations...,result.converged)
            println(io,join(values,'\t')); flush(io)
            path = joinpath(countsdir,"julia_$(Threads.nthreads())_$(mode)_$(rep).tsv")
            open(path,"w") do out
                println(out,"sequence\tcount")
                for j in eachindex(table.sequences)
                    println(out,table.sequences[j],'\t',table.counts[1,j])
                end
            end
            println("Julia ",Threads.nthreads()," threads ",mode," rep ",rep,": ",round(total;digits=3)," s")
            flush(stdout)
        end
    end
end
main()
