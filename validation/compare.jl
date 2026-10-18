using DADA2

root = joinpath(@__DIR__,"generated")
function read_errors(path)
    rows = [parse.(Float64,split(line,'\t')) for line in readlines(path)]
    ErrorModel(permutedims(hcat(rows...)))
end
results = DadaResult[]
for direction in ("forward","reverse")
    err = read_errors(joinpath(root,"errors_$direction.tsv"))
    result = dada(derepFastq(joinpath(root,"$direction.fastq.gz")),err)
    push!(results,result)
    expected = Dict{String,Int}()
    for row in readlines(joinpath(root,"asvs_$direction.tsv"))[2:end]
        s,c = split(row,'\t'); expected[s] = parse(Int,c)
    end
    actual = Dict(zip(result.sequences,result.counts))
    keys_union = union(keys(expected),keys(actual))
    shared = intersect(keys(expected),keys(actual))
    similarity = sum(min(get(expected,s,0),get(actual,s,0)) for s in keys_union) /
                 sum(max(get(expected,s,0),get(actual,s,0)) for s in keys_union)
    println(direction,": R_ASVs=",length(expected)," Julia_ASVs=",length(actual),
            " shared_ASVs=",length(shared)," abundance_weighted_Jaccard=",round(similarity;digits=6),
            " exact_counts=",expected==actual," R_reads=",sum(values(expected)),
            " Julia_reads=",sum(values(actual))," converged=",result.converged)
    @assert sum(values(actual)) + result.rejected == sum(result.derep.counts)
    @assert actual == expected "Inference differs from the pinned R fixture"
end
merged = mergePairs(results[1],results[2])
table = removeBimeraDenovo(makeSequenceTable([merged]))
expected = Dict{String,Int}()
for row in readlines(joinpath(root,"asvs_final.tsv"))[2:end]
    s,c = split(row,'\t'); expected[s] = parse(Int,c)
end
actual = Dict(zip(table.sequences,vec(table.counts)))
println("final: R_ASVs=",length(expected)," Julia_ASVs=",length(actual),
        " exact_counts=",actual==expected," R_reads=",sum(values(expected)),
        " Julia_reads=",sum(values(actual)))
@assert actual == expected "Final sequence table differs from the pinned R fixture"
