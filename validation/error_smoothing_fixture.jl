using DADA2

# Curved log rates, unequal exposures and unobserved tails: no DADA2 inference.
trans = zeros(Int,16,51)
for from in 1:4, q in 10:40
    total = 2_000_000 + 50_000*q + 100_000*from
    for to in 1:4
        from == to && continue
        rate = 10.0^(-2.5 - 0.04*q - 0.001*(q-25)^2 + 0.04*to)
        trans[4(from-1)+to,q+1] = round(Int,total*rate)
    end
    trans[4(from-1)+from,q+1] = total-sum(trans[4(from-1)+to,q+1] for to in 1:4 if to != from)
end
path = joinpath(@__DIR__,"..","test","fixtures","error_transitions.tsv")
open(path,"w") do io
    for row in eachrow(trans)
        println(io,join(row,'\t'))
    end
end
