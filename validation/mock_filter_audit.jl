using DADA2
root = joinpath(@__DIR__,"generated","public")
function paired_records(f,r)
    records = Dict{String,Tuple{FastqRecord,FastqRecord}}()
    fr = FastqReader(f); rr = FastqReader(r)
    try
        for (a,b) in zip(fr,rr)
            records[first(split(a.id))] = (a,b)
        end
    finally
        close(fr); close(rr)
    end
    records
end
j = paired_records(joinpath(root,"j_filteredF.fastq.gz"),joinpath(root,"j_filteredR.fastq.gz"))
r = paired_records(joinpath(root,"r_filteredF.fastq.gz"),joinpath(root,"r_filteredR.fastq.gz"))
shared = intersect(keys(j),keys(r))
same = count(id -> j[id][1].sequence == r[id][1].sequence && j[id][2].sequence == r[id][2].sequence &&
                  j[id][1].quality == r[id][1].quality && j[id][2].quality == r[id][2].quality,shared)
open(joinpath(@__DIR__,"mock_filter_audit.tsv"),"w") do io
    println(io,"id\tpresent_in\tforward_length\treverse_length")
    for (label,current,other) in (("Julia_only",j,r),("R_only",r,j))
        for id in sort!(collect(setdiff(keys(current),keys(other))))
            a,b = current[id]
            println(io,id,'\t',label,'\t',length(a.sequence),'\t',length(b.sequence))
        end
    end
end
println("Shared pairs: ",length(shared),"; identical sequence/quality pairs: ",same,
        "; Julia-only: ",length(setdiff(keys(j),keys(r))),"; R-only: ",length(setdiff(keys(r),keys(j))))
