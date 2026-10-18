using DADA2, Test

root = joinpath(@__DIR__,"generated","public")
f = derepFastq(joinpath(root,"r_filteredF.fastq.gz"))
r = derepFastq(joinpath(root,"r_filteredR.fastq.gz"))
references = readFasta(joinpath(root,"HMP_MOCK.v35.fasta")).sequences
open(joinpath(@__DIR__,"error_learning_mock.tsv"),"w") do io
    println(io,"method\tasvs\texact_matches\tunmatched\treads\tforward_converged\treverse_converged")
    for smoother in (:linear,:loess,:raw)
        ef = smoother == :raw ? learnErrors(joinpath(root,"r_filteredF.fastq.gz")) : learnErrors(f;smoother)
        er = smoother == :raw ? learnErrors(joinpath(root,"r_filteredR.fastq.gz")) : learnErrors(r;smoother)
        merged = mergePairs(dada(f,ef),dada(r,er))
        # Same expected-length selection as the frozen mock protocol.
        table = makeSequenceTable([merged])
        keep = findall(s -> 250 <= length(s) <= 256,table.sequences)
        table = SequenceTable(table.samples,table.sequences[keep],table.counts[:,keep])
        clean = removeBimeraDenovo(table;method=:pooled)
        matched = count(s -> any(ref -> occursin(s,ref) || occursin(reverseComplement(s),ref),references),clean.sequences)
        println(io,join((smoother,length(clean.sequences),matched,length(clean.sequences)-matched,
            sum(clean.counts),ef.converged,er.converged),'\t')); flush(io)
        @test matched == length(clean.sequences) == 20
        @test ef.converged && er.converged
    end
end
