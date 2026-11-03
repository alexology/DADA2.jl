# Compare the outputs of benchmark/pipeline_timing.R and .jl in `out`.
# Usage: julia --project=. benchmark/pipeline_compare.jl out
using CodecZlib

out = only(ARGS)
R = joinpath(out, "R"); J = joinpath(out, "Julia")
table(path) = Dict(String(a) => parse(Int, b) for (a, b) in (split(l, '\t') for l in readlines(path)[2:end]))
matrix(path) = permutedims(hcat([parse.(Float64, split(l, '\t')) for l in readlines(path)]...))
reads(path) = open(io -> readlines(GzipDecompressorStream(io)), path)
for side in ("F", "R")
    er = matrix(joinpath(R, "err_$side.tsv")); ej = matrix(joinpath(J, "err_$side.tsv"))
    difference = size(er) == size(ej) ? maximum(abs.(er .- ej) ./ er) : NaN
    println(side, ": filtered reads identical=",
            reads(joinpath(R, "$side.fastq.gz")) == reads(joinpath(J, "$side.fastq.gz")),
            ", error model max relative difference=", difference)
end
for name in ("asv_F", "asv_R", "merged", "table", "nochim")
    a = table(joinpath(R, "$name.tsv")); b = table(joinpath(J, "$name.tsv"))
    println(rpad(name, 8), " R=", length(a), " Julia=", length(b),
            " shared=", length(intersect(keys(a), keys(b))),
            " reads R/Julia=", sum(values(a); init=0), "/", sum(values(b); init=0), " identical=", a == b)
end
for side in ("R", "Julia")
    println("\n", side, ":\n", read(joinpath(out, side, "summary.txt"), String))
end
