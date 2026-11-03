# Paired-end pipeline in DADA2.jl, timed by step, for comparison with
# benchmark/pipeline_timing.R (same settings and outputs). The code is first
# compiled on the bundled fixture, so that compilation is not timed.
# Usage: julia --project=. --threads=4 benchmark/pipeline_timing.jl R1.fastq.gz R2.fastq.gz out \
#          [trimLeftF trimLeftR truncLenF truncLenR maxEE]
using DADA2

length(ARGS) >= 3 || error("Usage: pipeline_timing.jl R1 R2 out [trimF trimR truncF truncR maxEE]")
arg(i, default) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : default
trimF, trimR, truncF, truncR = Int(arg(4, 0)), Int(arg(5, 0)), Int(arg(6, 0)), Int(arg(7, 0))
maxee = arg(8, 2.0)

function write_table(path, sequences, counts)
    open(path, "w") do io
        println(io, "sequence\tabundance")
        foreach((s, c) -> println(io, s, '\t', c), sequences, counts)
    end
end

function pipeline(r1, r2, out, optsF, optsR; save=true)
    mkpath(out)
    f = joinpath(out, "F.fastq.gz"); r = joinpath(out, "R.fastq.gz")
    times = Pair{String,Float64}[]
    t = @elapsed stats = filterAndTrim(r1, r2, f, r; optionsF=optsF, optionsR=optsR)
    push!(times, "filter" => t)
    t = @elapsed (ef = learnErrors(f); er = learnErrors(r))
    push!(times, "learn" => t)
    t = @elapsed (df = derepFastq(f); dr = derepFastq(r))
    push!(times, "derep" => t)
    t = @elapsed (af = dada(df, ef); ar = dada(dr, er))
    push!(times, "dada" => t)
    t = @elapsed merged = mergePairs(af, ar)
    push!(times, "merge" => t)
    t = @elapsed (table = makeSequenceTable([merged]); nochim = removeBimeraDenovo(table))
    push!(times, "chimeras" => t)
    save || return
    for (name, e) in (("err_F.tsv", ef), ("err_R.tsv", er))
        open(joinpath(out, name), "w") do io
            foreach(row -> println(io, join(repr.(row), '\t')), eachrow(e.errors.rates))
        end
    end
    write_table(joinpath(out, "asv_F.tsv"), af.sequences, af.counts)
    write_table(joinpath(out, "asv_R.tsv"), ar.sequences, ar.counts)
    write_table(joinpath(out, "merged.tsv"), merged.sequences, merged.counts)
    write_table(joinpath(out, "table.tsv"), table.sequences, vec(table.counts))
    write_table(joinpath(out, "nochim.tsv"), nochim.sequences, vec(nochim.counts))
    open(joinpath(out, "summary.txt"), "w") do io
        println(io, "reads_in ", stats.reads_in, "\nreads_out ", stats.reads_out)
        println(io, "rounds_F ", ef.iterations, "\nrounds_R ", er.iterations)
        foreach(((k, v),) -> println(io, "time_", k, " ", round(v; digits=2)), times)
    end
    print(read(joinpath(out, "summary.txt"), String))
end

fixtures = joinpath(pkgdir(DADA2), "test", "fixtures")
mktempdir() do dir
    pipeline(joinpath(fixtures, "sam1F.fastq.gz"), joinpath(fixtures, "sam1R.fastq.gz"), dir,
             FilterOptions(truncLen=240, maxEE=2.0), FilterOptions(truncLen=160, maxEE=2.0); save=false)
end
println("Julia ", VERSION, ", threads=", Threads.nthreads())
pipeline(ARGS[1], ARGS[2], joinpath(ARGS[3], "Julia"),
         FilterOptions(trimLeft=trimF, truncLen=truncF, maxEE=maxee),
         FilterOptions(trimLeft=trimR, truncLen=truncR, maxEE=maxee))
