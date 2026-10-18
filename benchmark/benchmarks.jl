using DADA2, Random

function dense_alignment_score(a,b)
    scores = zeros(Int,length(a)+1,length(b)+1)
    scores[:,1] .= -8 .* (0:length(a)); scores[1,:] .= -8 .* (0:length(b))
    for i in 1:length(a), j in 1:length(b)
        scores[i+1,j+1] = max(scores[i,j]+(a[i]==b[j] ? 5 : -4),scores[i,j+1]-8,scores[i+1,j]-8)
    end
    scores[end,end]
end
function measure(f; repeats=5)
    f() # compile before timing
    times = Float64[]; bytes = Int[]
    for _ in 1:repeats
        GC.gc()
        result = @timed f()
        push!(times,result.time); push!(bytes,result.bytes)
    end
    minimum(times),minimum(bytes)
end
rng = MersenneTwister(42)
a = String(rand(rng,UInt8.(['A','C','G','T']),1000))
bytes = collect(codeunits(a)); bytes[500] = bytes[500] == 0x41 ? 0x43 : 0x41; b = String(bytes)
@assert alignSequences(a,b).score == dense_alignment_score(a,b)
rows = [("alignment_banded",measure(() -> alignSequences(a,b))),
        ("alignment_dense_score_only",measure(() -> dense_alignment_score(a,b)))]
seqs = [String(rand(rng,UInt8.(['A','C','G','T']),250)) for _ in 1:100]
d = Derep(seqs,fill(10,100),[fill(35.0,250) for _ in seqs])
push!(rows,("inference_serial",measure(() -> dada(d;threaded=false))))
push!(rows,("inference_threaded",measure(() -> dada(d;threaded=true))))
mktempdir() do dir
    reference = joinpath(dir,"taxonomy.fa")
    writeFasta(reference,seqs;ids=["Bacteria;Taxon$i;" for i in eachindex(seqs)])
    model = trainTaxonomy(reference)
    push!(rows,("taxonomy_train_and_classify",measure(() -> assignTaxonomy(seqs[1:10],reference;nBoot=10))))
    push!(rows,("taxonomy_reuse_index",measure(() -> assignTaxonomy(seqs[1:10],model;nBoot=10))))
end
println("Julia ",VERSION,"; threads=",Threads.nthreads())
open(joinpath(@__DIR__,"results.tsv"),"w") do io
    println(io,"case\tseconds\tallocated_bytes")
    for (name,(seconds,bytes)) in rows
        println(io,name,'\t',seconds,'\t',bytes)
        println(name,": ",round(seconds;digits=6)," s, ",bytes," allocated bytes")
    end
end
