using DADA2, Test

# Optional integration check against the pinned upstream FASTQ fixtures.
root = joinpath(@__DIR__,"fixtures")
isdir(root) || error("Missing bundled FASTQ fixtures")
mktempdir() do dir
    f = joinpath(dir,"f.fastq.gz"); r = joinpath(dir,"r.fastq.gz")
    stats = filterAndTrim(joinpath(root,"sam1F.fastq.gz"),joinpath(root,"sam1R.fastq.gz"),f,r;
        optionsF=FilterOptions(truncLen=240,maxEE=2.0),
        optionsR=FilterOptions(truncLen=160,maxEE=2.0))
    fd = derepFastq(f); rd = derepFastq(r)
    ef = learnErrors(fd;maxIter=10); er = learnErrors(rd;maxIter=10)
    fr = dada(fd,ef); rr = dada(rd,er)
    merged = mergePairs(fr,rr)
    table = removeBimeraDenovo(makeSequenceTable([merged]))
    @test sum(fr.counts) + fr.rejected == sum(rr.counts) + rr.rejected == stats.reads_out
    @test merged.accepted + merged.rejected == stats.reads_out
    @test !isempty(table.sequences)
    @test sum(table.counts) <= merged.accepted
    println("reads_in=",stats.reads_in," reads_out=",stats.reads_out,
            " forward_ASVs=",length(fr.sequences)," reverse_ASVs=",length(rr.sequences),
            " merged_reads=",merged.accepted," merged_ASVs=",length(merged.sequences),
            " nonchimeric_ASVs=",length(table.sequences),
            " error_converged=",(ef.converged,er.converged),
            " inference_converged=",(fr.converged,rr.converged))
end
