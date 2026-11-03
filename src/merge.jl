# Port of R dada2 `mergePairs` (paired.R) with `C_eval_pair` and
# `C_pair_consensus` (evaluate.cpp). Alignment rows use 0 for gaps.

# Row access in which positions past either end behave like C's '\0'.
@inline row_at(v, i) = 0 <= i < length(v) ? v[i+1] : 0xff

"Upstream `C_eval_pair`: matches, mismatches and indels inside the overlap."
function eval_pair(al0::Vector{UInt8}, al1::Vector{UInt8})
    n = length(al0)
    g1 = g2 = true; start = -1
    while true
        start += 1
        g1 = g1 && row_at(al0, start) == 0x00; g2 = g2 && row_at(al1, start) == 0x00
        ((g1 || g2) && start < n) || break
    end
    g1 = g2 = true; stop = n
    while true
        stop -= 1
        g1 = g1 && row_at(al0, stop) == 0x00; g2 = g2 && row_at(al1, stop) == 0x00
        ((g1 || g2) && stop >= start) || break
    end
    match = mismatch = indel = 0
    for i in start:stop
        a = row_at(al0, i); b = row_at(al1, i)
        if a == 0x00 || b == 0x00
            indel += 1
        elseif a == b
            match += 1
        else
            mismatch += 1
        end
    end
    match, mismatch, indel
end

"Upstream `C_pair_consensus`: overlap consensus, mismatches from the preferred read."
function pair_consensus(al0::Vector{UInt8}, al1::Vector{UInt8}, prefer::Int, trim_overhang::Bool)
    n = length(al0)
    out = Vector{UInt8}(undef, n)
    for i in 1:n
        a = al0[i]; b = al1[i]
        out[i] = a == b || b == 0x00 ? a : a == 0x00 ? b : prefer == 1 ? a : prefer == 2 ? b : 0x05
    end
    if trim_overhang
        for i in 1:n
            al0[i] != 0x00 && break
            out[i] = 0x00
        end
        for i in n:-1:1
            al1[i] != 0x00 && break
            out[i] = 0x00
        end
    end
    String([b == 0x05 ? 'N' : "ACGT"[b] for b in out if b != 0x00])
end

"""
    mergePairs(forward, reverse; minOverlap=12, maxMismatch=0, justConcatenate=false, trimOverhang=false)

Port of R `mergePairs`. Each pair of forward/reverse ASVs observed in the reads
is aligned (unbanded ends-free Needleman–Wunsch of the forward ASV and the
reverse complement of the reverse ASV; match 1, mismatch and gap -64, or -8
when `maxMismatch > 0`) and accepted with at least `minOverlap` matches and at
most `maxMismatch` mismatches plus indels in the overlap. Mismatches take the
base of the read whose ASV has more reads without substitutions (`n0`).
`justConcatenate` joins the reads with ten `N`; `trimOverhang` removes the
parts that extend past the start of the other read. Requires the read maps
(`keepmap=true`).
"""
function mergePairs(forward::DadaResult, reverse::DadaResult; minOverlap=12, maxMismatch=0,
                    justConcatenate=false, trimOverhang=false)
    minOverlap >= 0 && maxMismatch >= 0 || throw(ArgumentError("Invalid merge options"))
    fm = forward.derep.map; rm = reverse.derep.map
    isempty(fm) && throw(ArgumentError("mergePairs requires keepmap=true"))
    length(fm) == length(rm) || throw(ArgumentError("Paired read counts differ"))
    order = Tuple{Int,Int}[]; abundance = Dict{Tuple{Int,Int},Int}()
    for i in eachindex(fm)
        key = (forward.assignments[fm[i]], reverse.assignments[rm[i]])
        (key[1] == 0 || key[2] == 0) && continue
        haskey(abundance, key) || push!(order, key)
        abundance[key] = get(abundance, key, 0) + 1
    end
    match, mismatch, gap = maxMismatch == 0 ? (1, -64, -64) : (1, -8, -8)
    rows = map(order) do (f, r)
        fseq = forward.sequences[f]; rseq = reverseComplement(reverse.sequences[r])
        justConcatenate && return (fseq * "NNNNNNNNNN" * rseq, 0, 0, 0, missing, true)
        al = nwalign_endsfree(encode_dna(fseq), encode_dna(rseq), match, mismatch, gap, -1)
        nmatch, nmismatch, nindel = eval_pair(al.al0, al.al1)
        prefer = 1 + (reverse.clustering.n0[r] > forward.clustering.n0[f])
        accept = nmatch >= minOverlap && nmismatch + nindel <= maxMismatch
        sequence = accept ? pair_consensus(al.al0, al.al1, prefer, trimOverhang) : ""
        (sequence, nmatch, nmismatch, nindel, prefer, accept)
    end
    ranked = sortperm(eachindex(order); by=k -> -abundance[order[k]])  # stable, as R's order()
    column(i) = [rows[k][i] for k in ranked]
    pairs = MergeTable((column(1), [abundance[order[k]] for k in ranked], [order[k][1] for k in ranked],
                        [order[k][2] for k in ranked], column(2), column(3), column(4),
                        Union{Missing,Int}[rows[k][5] for k in ranked], column(6)))
    # Identical merged sequences are summed, as R's getUniques.
    merged = Dict{String,Int}()
    for k in eachindex(pairs.sequence)
        pairs.accept[k] || continue
        merged[pairs.sequence[k]] = get(merged, pairs.sequence[k], 0) + pairs.abundance[k]
    end
    sequences = sort!(collect(keys(merged)); by=s -> (-merged[s], s))
    accepted = sum(values(merged); init=0)
    MergeResult(sequences, [merged[s] for s in sequences], accepted, length(fm) - accepted, pairs)
end
