using DADA2

root = joinpath(@__DIR__,"generated","public")
refs = readFasta(joinpath(root,"HMP_MOCK.v35.fasta"))
function readerr(path)
    rows = [parse.(Float64,split(line,'\t')) for line in readlines(path)]
    ErrorModel(permutedims(hcat(rows...)))
end
function readcounts(path)
    pairs = [split(line,'\t') for line in readlines(path)[2:end]]
    Dict(s=>parse(Int,n) for (s,n) in pairs)
end
results = Dict{String,Dict{String,Int}}("r_learned"=>readcounts(joinpath(root,"r_mock.tsv")))
fpath = joinpath(root,"j_filteredF.fastq.gz"); rpath = joinpath(root,"j_filteredR.fastq.gz")
if !isfile(fpath) && !isfile(rpath)
    stats = filterAndTrim(joinpath(root,"Mock_S280_L001_R1_001.fastq"),joinpath(root,"Mock_S280_L001_R2_001.fastq"),fpath,rpath;
        optionsF=FilterOptions(truncLen=240,maxEE=2.0),optionsR=FilterOptions(truncLen=160,maxEE=2.0))
    println("Julia filter: ",stats.reads_in," → ",stats.reads_out)
end
for mode in ("r_errors","learned")
    f = derepFastq(mode=="r_errors" ? joinpath(root,"r_filteredF.fastq.gz") : fpath)
    r = derepFastq(mode=="r_errors" ? joinpath(root,"r_filteredR.fastq.gz") : rpath)
    if mode == "learned"
        ef = learnErrors(f;maxIter=10); er = learnErrors(r;maxIter=10)
        println("Julia mock error convergence: ",(ef.converged,er.converged))
        errf = ef.errors; errr = er.errors
    else
        errf = readerr(joinpath(root,"mock_errF.tsv")); errr = readerr(joinpath(root,"mock_errR.tsv"))
    end
    df = dada(f,errf); dr = dada(r,errr)
    tab = makeSequenceTable([mergePairs(df,dr)])
    keep = [250 <= ncodeunits(s) <= 256 for s in tab.sequences]
    tab = SequenceTable(tab.samples,tab.sequences[keep],tab.counts[:,keep])
    tab = removeBimeraDenovo(tab;method=:pooled)
    actual = Dict(zip(tab.sequences,vec(tab.counts)))
    results["julia_$mode"] = actual
    open(joinpath(root,"julia_mock_$mode.tsv"),"w") do io
        println(io,"sequence\tcount")
        for s in sort!(collect(keys(actual)))
            println(io,s,'\t',actual[s])
        end
    end
end
open(joinpath(@__DIR__,"mock_results.tsv"),"w") do out
    println(out,"method\tasvs\texact_reference_matches\tunmatched_asvs\texact_match_precision\tunmatched_read_fraction\taccepted_reads\treference_records_covered")
    for method in sort!(collect(keys(results)))
        counts = results[method]
        matched = [s for s in keys(counts) if any(r -> occursin(s,r) || occursin(reverseComplement(s),r),refs.sequences)]
        total = sum(values(counts))
        unmatched = sum(counts[s] for s in setdiff(keys(counts),matched);init=0)
        covered = count(r -> any(s -> occursin(s,r) || occursin(reverseComplement(s),r),keys(counts)),refs.sequences)
        row = (method,length(counts),length(matched),length(counts)-length(matched),
               length(matched)/length(counts),unmatched/total,total,covered)
        println(out,join(row,'\t')); println(join(row,'\t'))
    end
end
println("Reference records: ",length(refs.sequences))
