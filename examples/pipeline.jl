using DADA2

# Run: julia --project=. --threads=auto examples/pipeline.jl reads.fastq output_directory
length(ARGS) == 2 || error("Usage: pipeline.jl reads.fastq output_directory")
input, output = ARGS
mkpath(output)
filtered = joinpath(output,"filtered.fastq.gz")
stats = filterAndTrim(input,filtered;options=FilterOptions(maxEE=2.0,minLen=20))
learned = learnErrors(filtered)
reads = derepFastq(filtered;keepmap=false)
learned.converged || @warn "Error learning did not converge" iterations=learned.iterations
result = dada(reads,learned)
result.converged || @warn "Variant inference did not converge"
table = removeBimeraDenovo(makeSequenceTable([result]))
writeFasta(joinpath(output,"asvs.fasta"),table.sequences)
open(joinpath(output,"counts.tsv"),"w") do io
    println(io,"ASV\tcount")
    for j in eachindex(table.sequences)
        println(io,"ASV",j,'\t',table.counts[1,j])
    end
end
println("Reads: ",stats.reads_in," → ",stats.reads_out,"; ASVs: ",length(table.sequences))
