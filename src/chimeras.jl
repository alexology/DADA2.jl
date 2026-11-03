# Port of R dada2 chimera detection (chimeras.R, chimera.cpp): alignments are
# the upstream vectorized ends-free kernel, banded by `maxShift`.

"Upstream `get_lr`: query bases explained from the left and from the right by a parent."
function get_lr(al0::Vector{UInt8}, al1::Vector{UInt8}, allow_one_off::Bool, max_shift::Int)
    len = length(al0)
    pos = 0; left = 0; left_oo = right_oo = 0
    while row_at(al0, pos) == 0x00 && pos < len
        pos += 1
    end
    while row_at(al1, pos) == 0x00 && pos < max_shift
        pos += 1; left += 1
    end
    while pos < len && al0[pos+1] == al1[pos+1]
        pos += 1; left += 1
    end
    if allow_one_off
        left_oo = left; pos += 1
        pos < len && al0[pos+1] != 0x00 && (left_oo += 1)
        while pos < len && al0[pos+1] == al1[pos+1]
            pos += 1; left_oo += 1
        end
    end
    pos = len - 1; right = 0
    while row_at(al0, pos) == 0x00 && pos >= 0
        pos -= 1
    end
    # Upstream compares `pos > len - max_shift` as size_t.
    while row_at(al1, pos) == 0x00 && (pos % UInt) > ((len - max_shift) % UInt)
        pos -= 1; right += 1
    end
    while pos >= 0 && al0[pos+1] == al1[pos+1]
        pos -= 1; right += 1
    end
    if allow_one_off
        right_oo = right; pos -= 1
        pos >= 0 && al0[pos+1] != 0x00 && (right_oo += 1)
        while pos >= 0 && al0[pos+1] == al1[pos+1]
            pos -= 1; right_oo += 1
        end
    end
    left, right, left_oo, right_oo
end

"Upstream `get_ham_endsfree`: differences inside the overlap."
function ham_endsfree(al0::Vector{UInt8}, al1::Vector{UInt8})
    i = 0
    g1 = row_at(al0, i) == 0x00; g2 = row_at(al1, i) == 0x00
    while g1 || g2
        i += 1
        g1 = g1 && row_at(al0, i) == 0x00; g2 = g2 && row_at(al1, i) == 0x00
    end
    j = length(al1) - 1
    g1 = row_at(al0, j) == 0x00; g2 = row_at(al1, j) == 0x00
    while g1 || g2
        j -= 1
        g1 = g1 && row_at(al0, j) == 0x00; g2 = g2 && row_at(al1, j) == 0x00
    end
    count(p -> row_at(al0, p) != row_at(al1, p), i:j)
end

struct ChimeraOptions
    allow_one_off::Bool
    min_one_off_par_dist::Int
    match::Int
    mismatch::Int
    gap::Int
    max_shift::Int
end
ChimeraOptions(; allowOneOff=false, minOneOffParentDistance=4, maxShift=16, match=5, mismatch=-4,
               gap=-8) =
    ChimeraOptions(allowOneOff, minOneOffParentDistance, match, mismatch, -abs(gap), maxShift)

function parent_overlap(query::Vector{UInt8}, parent::Vector{UInt8}, o::ChimeraOptions, scratch)
    al = nwalign_vectorized(query, parent, o.match, o.mismatch, o.gap, 0, o.max_shift; scratch)
    left, right, left_oo, right_oo = get_lr(al.al0, al.al1, o.allow_one_off, o.max_shift)
    allowed = o.allow_one_off && ham_endsfree(al.al0, al.al1) >= o.min_one_off_par_dist
    left, right, left_oo, right_oo, allowed
end

"Upstream `C_is_bimera`."
function is_bimera(query::Vector{UInt8}, parents, o::ChimeraOptions, scratch)
    n = length(query)
    max_left = max_right = 0
    oo_left = oo_right = oo_left_oo = oo_right_oo = 0
    for parent in parents
        left, right, left_oo, right_oo, allowed = parent_overlap(query, parent, o, scratch)
        left + right >= n && continue  # identical, pure shift or internal indel
        max_left = max(max_left, left); max_right = max(max_right, right)
        if allowed
            oo_left = max(oo_left, left); oo_right = max(oo_right, right)
            oo_left_oo = max(oo_left_oo, left_oo); oo_right_oo = max(oo_right_oo, right_oo)
        end
        max_left + max_right >= n && return true
        o.allow_one_off && (oo_left + oo_right_oo >= n || oo_left_oo + oo_right >= n) && return true
    end
    false
end

"""
    isBimera(sequence, parents; allowOneOff=false, minOneOffParentDistance=4, maxShift=16)

Port of R `isBimera`: true if the sequence is an exact (or, with
`allowOneOff`, one-mismatch) join of a left part of one parent and a right
part of another, aligned ends-free with shifts up to `maxShift`; parents that
explain the whole sequence alone (identical, shifted, internal indel) are
ignored.
"""
isBimera(sequence::DNAInput, parents::AbstractVector{<:DNAInput}; kwargs...) =
    is_bimera(encode_dna(validated_dna(sequence)), [encode_dna(validated_dna(p)) for p in parents],
              ChimeraOptions(; kwargs...), AlignScratch())

"""
    isBimeraDenovo(sequences, abundances; minFoldParentOverAbundance=2, minParentAbundance=8, ...)
    isBimeraDenovo(table; ...)

Port of R `isBimeraDenovo`: each sequence is tested against those more than
`minFoldParentOverAbundance` times and more than `minParentAbundance` reads
abundant. A `SequenceTable` is tested on column totals.
"""
function isBimeraDenovo(sequences::AbstractVector{<:DNAInput}, abundances::AbstractVector{<:Integer};
                        minFoldParentOverAbundance=2, minParentAbundance=8, threaded=true, kwargs...)
    length(sequences) == length(abundances) || throw(ArgumentError("One abundance per sequence is required"))
    o = ChimeraOptions(; kwargs...)
    seqs = [encode_dna(validated_dna(s)) for s in sequences]
    flags = fill(false, length(seqs))  # not a BitVector: written by several tasks
    foreach_scratch(length(seqs), threaded) do scratch, j
        a = abundances[j]
        parents = [seqs[k] for k in eachindex(seqs)
                   if abundances[k] > minFoldParentOverAbundance*a && abundances[k] > minParentAbundance]
        length(parents) >= 2 && (flags[j] = is_bimera(seqs[j], parents, o, scratch))
    end
    flags
end
isBimeraDenovo(table::SequenceTable; kwargs...) =
    isBimeraDenovo(table.sequences, vec(sum(table.counts; dims=1)); kwargs...)

"""
    isBimeraDenovoTable(table; minSampleFraction=0.9, ignoreNNegatives=1,
                        minFoldParentOverAbundance=1.5, minParentAbundance=2, ...)

Port of R `isBimeraDenovoTable` (consensus): a sequence is tested in every
sample containing it, with parents more than `minFoldParentOverAbundance`
times and at least `minParentAbundance` reads abundant in that sample. It is a
chimera if flagged in all those samples, or in at least `minSampleFraction`
of them after discounting `ignoreNNegatives`.
"""
function isBimeraDenovoTable(table::SequenceTable; minSampleFraction=0.9, ignoreNNegatives=1,
                             minFoldParentOverAbundance=1.5, minParentAbundance=2, threaded=true,
                             kwargs...)
    o = ChimeraOptions(; kwargs...)
    seqs = [encode_dna(validated_dna(s)) for s in table.sequences]
    counts = table.counts; ncol = length(seqs)
    nflag = zeros(Int, ncol); nsam = zeros(Int, ncol)
    foreach_scratch(ncol, threaded) do scratch, j
        # Alignments with each parent are computed once per sequence.
        overlaps = Vector{Union{Nothing,Tuple{Int,Int,Int,Int,Bool}}}(nothing, ncol)
        n = length(seqs[j])
        for i in axes(counts, 1)
            counts[i,j] <= 0 && continue
            nsam[j] += 1
            max_left = max_right = 0
            oo_left = oo_right = oo_left_oo = oo_right_oo = 0
            for k in 1:ncol
                counts[i,k] > minFoldParentOverAbundance*counts[i,j] &&
                    counts[i,k] >= minParentAbundance || continue
                if overlaps[k] === nothing
                    left, right, left_oo, right_oo, allowed = parent_overlap(seqs[j], seqs[k], o, scratch)
                    # Identical, pure-shift and internal-indel parents count as zero.
                    overlaps[k] = left + right < n ? (left, right, left_oo, right_oo, allowed) :
                                                     (0, 0, 0, 0, allowed)
                end
                left, right, left_oo, right_oo, allowed = overlaps[k]
                max_left = max(max_left, left); max_right = max(max_right, right)
                if o.allow_one_off && allowed
                    oo_left = max(oo_left, left); oo_right = max(oo_right, right)
                    oo_left_oo = max(oo_left_oo, left_oo); oo_right_oo = max(oo_right_oo, right_oo)
                end
            end
            if max_left + max_right >= n
                nflag[j] += 1
            elseif o.allow_one_off && (oo_left + oo_right_oo >= n || oo_left_oo + oo_right >= n)
                nflag[j] += 1
            end
        end
    end
    [nflag[j] >= nsam[j] || (nflag[j] > 0 && nflag[j] >= (nsam[j] - ignoreNNegatives)*minSampleFraction)
     for j in 1:ncol]
end

"""
    removeBimeraDenovo(table; method=:consensus, kwargs...)

Port of R `removeBimeraDenovo` for sequence tables: `:consensus`
(`isBimeraDenovoTable`), `:pooled` (`isBimeraDenovo` on column totals) or
`:per_sample` (chimeras zeroed sample by sample, empty columns dropped).
Other keywords go to the chosen test.
"""
function removeBimeraDenovo(table::SequenceTable; method=:consensus, kwargs...)
    if method == :consensus || method == :pooled
        bim = method == :consensus ? isBimeraDenovoTable(table; kwargs...) : isBimeraDenovo(table; kwargs...)
        keep = .!bim
        return SequenceTable(table.samples, table.sequences[keep], table.counts[:, keep])
    end
    method == :per_sample || throw(ArgumentError("method must be :consensus, :pooled or :per_sample"))
    counts = copy(table.counts)
    for i in axes(counts, 1)
        counts[i, isBimeraDenovo(table.sequences, counts[i, :]; kwargs...)] .= 0
    end
    keep = vec(sum(counts; dims=1)) .> 0
    SequenceTable(table.samples, table.sequences[keep], counts[:, keep])
end
