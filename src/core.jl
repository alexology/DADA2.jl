# Port of the DADA2 1.36.0 C++ inference core: Rmain.cpp, cluster.cpp,
# containers.cpp, pval.cpp, kmers.cpp, nwalign_endsfree.cpp and
# nwalign_vectorized.cpp. Iteration order, tie-breaking and the order of
# floating-point operations follow upstream, so that partitions are identical
# to R for the same error model. Indices are 1-based; sequences are encoded as
# A=1, C=2, G=3, T=4 and gaps as 0.

const KMER_SIZE = 5
const MAX_SHUFFLE = 10
const TAIL_APPROX_CUTOFF = 1e-7

encode_base(b::UInt8) = baseindex(b) % UInt8
encode_dna(s::AbstractString) = UInt8[encode_base(b) for b in codeunits(s)]

# ---------------------------------------------------------------------------
# Alignment (nwalign_endsfree.cpp, nwalign_vectorized.cpp)

"Alignment rows of the reference and query; zero marks a gap."
struct PairAlignment
    al0::Vector{UInt8}
    al1::Vector{UInt8}
end

"Reusable buffers for one task; an alignment returned from them is valid until the next call."
struct AlignScratch
    d::Vector{Int16}
    d32::Vector{Int32}
    p::Vector{UInt8}
    reversed::Vector{UInt8}  # first sequence reversed: anti-diagonals read it forwards
    al0::Vector{UInt8}
    al1::Vector{UInt8}
end
AlignScratch() = AlignScratch(Int16[], Int32[], UInt8[], UInt8[], UInt8[], UInt8[])
score_buffer(s::AlignScratch, ::Type{Int16}) = s.d
score_buffer(s::AlignScratch, ::Type{Int32}) = s.d32

function traceback_full(p, s1, s2, len1, len2)
    ncol = len2 + 1
    al0 = UInt8[]; al1 = UInt8[]
    i = len1; j = len2
    @inbounds while i > 0 || j > 0
        move = p[i*ncol+j+1]
        if move == 1
            i -= 1; j -= 1; push!(al0, s1[i+1]); push!(al1, s2[j+1])
        elseif move == 2
            j -= 1; push!(al0, 0x00); push!(al1, s2[j+1])
        elseif move == 3
            i -= 1; push!(al0, s1[i+1]); push!(al1, 0x00)
        else
            error("N-W Align out of range.")
        end
    end
    PairAlignment(reverse!(al0), reverse!(al1))
end

# Homopolymer flags: runs of at least three identical bases.
function homopolymer_flags(s)
    flags = falses(length(s))
    i = 1
    for j in eachindex(s)
        if j == length(s) || s[j] != s[j+1]
            j - i >= 2 && (flags[i:j] .= true)
            i = j + 1
        end
    end
    flags
end

"""
Banded ends-free Needleman–Wunsch on the full matrix (`nwalign_endsfree` and,
with `homo_gap != gap`, `nwalign_endsfree_homo`). `endsfree=false` gives the
global variant (`nwalign`). Ties prefer up, then left, then the diagonal.
"""
function nwalign_endsfree(s1, s2, match::Integer, mismatch::Integer, gap::Integer, band::Integer;
                          homo_gap::Integer=gap, endsfree::Bool=true)
    len1 = length(s1); len2 = length(s2)
    ncol = len2 + 1
    d = zeros(Int, (len1+1)*ncol); p = zeros(Int8, (len1+1)*ncol)
    homo = homo_gap != gap
    h1 = homo ? homopolymer_flags(s1) : falses(0)
    h2 = homo ? homopolymer_flags(s2) : falses(0)
    @inbounds for i in 0:len1
        d[i*ncol+1] = endsfree ? 0 : i*gap; p[i*ncol+1] = 3
    end
    @inbounds for j in 0:len2
        d[j+1] = endsfree ? 0 : j*gap; p[j+1] = 2
    end
    lband, rband = len2 > len1 ? (band, band+len2-len1) :
                   len1 > len2 ? (band+len1-len2, band) : (band, band)
    if band >= 0 && (band < len1 || band < len2)
        @inbounds for i in 0:len1
            i-lband-1 >= 0 && (d[i*ncol+i-lband-1+1] = -9999)
            i+rband+1 <= len2 && (d[i*ncol+i+rband+1+1] = -9999)
        end
    end
    @inbounds for i in 1:len1
        l, r = band >= 0 ? (max(1, i-lband), min(len2, i+rband)) : (1, len2)
        for j in l:r
            left = d[i*ncol+j] + (endsfree && i == len1 ? 0 : homo && h2[j] ? homo_gap : gap)
            up = d[(i-1)*ncol+j+1] + (endsfree && j == len2 ? 0 : homo && h1[i] ? homo_gap : gap)
            diag = d[(i-1)*ncol+j] + (s1[i] == s2[j] ? match : mismatch)
            if up >= diag && up >= left
                d[i*ncol+j+1] = up; p[i*ncol+j+1] = 3
            elseif left >= diag
                d[i*ncol+j+1] = left; p[i*ncol+j+1] = 2
            else
                d[i*ncol+j+1] = diag; p[i*ncol+j+1] = 1
            end
        end
    end
    traceback_full(p, s1, s2, len1, len2)
end

"Gapless alignment: positions paired from the start, the longer tail gapped."
function nwalign_gapless(s1, s2)
    n = max(length(s1), length(s2))
    PairAlignment(UInt8[i <= length(s1) ? s1[i] : 0x00 for i in 1:n],
                  UInt8[i <= length(s2) ? s2[i] : 0x00 for i in 1:n])
end

"""
Anti-diagonal banded aligner of upstream `nwalign_vectorized2`, translated
statement by statement. The shorter sequence is aligned as the first one;
`end_gap=0` makes the alignment ends-free.

Upstream scores in Int16, which overflows for long sequences (from about
6,500 identical bases with match 5; DD5): the traceback then loses pairs.
Int16 is kept whenever every score is bounded well inside its range, which
gives upstream's results; otherwise the same algorithm runs in Int32. The
two agree whenever Int16 does not overflow.
"""
function nwalign_vectorized(s1, s2, match::Integer, mismatch::Integer, gap::Integer, end_gap::Integer,
                            band::Integer; scratch::AlignScratch=AlignScratch(), compact::Val=Val(false))
    shorter = min(length(s1), length(s2)); total = length(s1) + length(s2)
    # In-band cells lie between the diagonal path from a free start and the
    # best all-match path; global end gaps add at most end_gap per base.
    bound = max(abs(match), abs(mismatch)) * shorter + (end_gap < 0 ? abs(end_gap) * total : 0) +
            4 * (abs(match) + abs(mismatch) + abs(gap))
    # Two static calls rather than a type chosen at run time (no dynamic dispatch).
    bound <= 32000 ?
        nwalign_vectorized(s1, s2, Int16(match), Int16(mismatch), Int16(gap), Int16(end_gap), Int(band), scratch, compact) :
        nwalign_vectorized(s1, s2, Int32(match), Int32(mismatch), Int32(gap), Int32(end_gap), Int(band), scratch, compact)
end

# With `Val(true)`, the traceback records only the substitutions (CompactSub),
# as needed while comparing with a centre; otherwise both alignment rows.
function nwalign_vectorized(s1::Vector{UInt8}, s2::Vector{UInt8}, match::T, mismatch::T,
                            gap::T, end_gap::T, band::Int, scratch::AlignScratch,
                            ::Val{compact}=Val(false)) where {T<:Union{Int16,Int32},compact}
    len1 = length(s1); len2 = length(s2)
    swap = false
    if len1 > len2
        s1, s2 = s2, s1; len1, len2 = len2, len1; swap = true
    end
    band < 0 && (band = len2)
    start_col = 1 + (1 + min(band, len1)) ÷ 2
    ncol = 2 + start_col + min(len2-len1+band, len2) ÷ 2
    nrow = len1 + len2 + 1
    fill_val = (Int(typemin(T)) - min(min(Int(mismatch), Int(gap)), min(Int(match), 0))) % T
    dbuffer = score_buffer(scratch, T)
    # As upstream, only the boundary columns are initialised: every other cell
    # is written before it is read (checked with randomised buffers against
    # the R oracles), so the buffers are reused without clearing.
    d = resize!(dbuffer, ncol*nrow)
    p = resize!(scratch.p, ncol*nrow)
    r1 = resize!(scratch.reversed, len1)
    @inbounds for i in 1:len1
        r1[i] = s1[len1-i+1]
    end
    at(row, col) = row*ncol + col + 1
    @inbounds for row in 0:nrow-1
        d[at(row,0)] = fill_val; d[at(row,1)] = fill_val
        d[at(row,ncol-2)] = fill_val; d[at(row,ncol-1)] = fill_val
    end
    d[at(0,start_col)] = zero(T)
    # Ends-free "left column" and "top row".
    row = 1; col = start_col-1; d_free = end_gap
    while row < 1 + min(band, len1)
        d[at(row,col)] = d_free; p[at(row,col)] = 3
        row % 2 == 0 && (col -= 1)
        row += 1; d_free += end_gap
    end
    row = 1; col = start_col; d_free = end_gap
    while row < 1 + min(band+len2-len1, len2)
        d[at(row,col)] = d_free; p[at(row,col)] = 2
        row % 2 == 1 && (col += 1)
        row += 1; d_free += end_gap
    end
    row = 2; col_min = start_col; col_max = start_col
    i_max = 0; j_min = 0; even = true
    recalc_left = false; recalc_right = false
    top = min(band, len1); right = min(band+len2-len1, len2)
    @inbounds while row <= len1 + len2
        width = col_max - col_min + 1
        previous = at(row-2, col_min)
        # s1[i_max-k+1] == r1[len1-i_max+k]: both sequences are read forwards.
        rbase = len1 - i_max; sbase = j_min + 1
        left_base = at(row-1, col_min - even)
        up_base = at(row-1, col_min + 1 - even)
        # Reads from rows row-2 and row-1, writes to row: no dependencies.
        out = at(row, col_min)
        if swap  # Precedence left > up > diag
            @simd ivdep for k in 0:width-1
                diag = d[previous+k] + ifelse(r1[rbase+k] == s2[sbase+k], match, mismatch)
                left = d[left_base+k] + gap; up = d[up_base+k] + gap
                first = left >= up
                entry = ifelse(first, left, up); pentry = ifelse(first, 0x02, 0x03)
                keep = entry >= diag
                d[out+k] = ifelse(keep, entry, diag); p[out+k] = ifelse(keep, pentry, 0x01)
            end
        else     # Precedence up > left > diag
            @simd ivdep for k in 0:width-1
                diag = d[previous+k] + ifelse(r1[rbase+k] == s2[sbase+k], match, mismatch)
                left = d[left_base+k] + gap; up = d[up_base+k] + gap
                first = up >= left
                entry = ifelse(first, up, left); pentry = ifelse(first, 0x03, 0x02)
                keep = entry >= diag
                d[out+k] = ifelse(keep, entry, diag); p[out+k] = ifelse(keep, pentry, 0x01)
            end
        end
        if row == top
            col_min -= 1; i_max += 1; j_min -= 1
        end
        row == right && (col_max += 1)
        if end_gap > gap
            if recalc_left
                d_free = d[left_base] + end_gap
                cell = at(row,col_min)
                if d_free > d[cell]
                    d[cell] = d_free; p[cell] = 2
                elseif !swap && d_free == d[cell] && p[cell] == 1
                    p[cell] = 2
                elseif swap && d_free == d[cell] && p[cell] != 2
                    p[cell] = 2
                end
            end
            i_max == len1 - 1 && (recalc_left = true)
            if recalc_right
                d_free = d[up_base + col_max - col_min] + end_gap
                cell = at(row,col_max)
                if d_free > d[cell]
                    d[cell] = d_free; p[cell] = 3
                elseif !swap && d_free == d[cell] && p[cell] != 3
                    p[cell] = 3
                elseif swap && d_free == d[cell] && p[cell] == 1
                    p[cell] = 3
                end
            end
            (row+1) ÷ 2 + col_max - start_col == len2 && (recalc_right = true)
        end
        if row < band && row < len1
            even && (col_min -= 1)
            i_max += 1
        elseif i_max < len1 - 1
            if band % 2 == 0
                even ? (j_min += 1) : (i_max += 1)
            elseif even
                col_min -= 1; i_max += 1
            else
                col_min += 1; j_min += 1
            end
        else
            even || (col_min += 1)
            j_min += 1
        end
        if row < right
            even || (col_max += 1)
        elseif (row+1) ÷ 2 + col_max - start_col < len2
            if (band+len2-len1) % 2 == 0
                even ? (col_max -= 1) : (col_max += 1)
            end
        else
            even && (col_max -= 1)
        end
        row += 1; even = !even
    end
    compact && return compact_traceback(p, s1, s2, len1, len2, ncol, start_col, swap)
    al0 = empty!(scratch.al0); al1 = empty!(scratch.al1)
    i = len1; j = len2
    @inbounds while i > 0 || j > 0
        offset = 2*start_col + j - i
        offset >= 0 || error("N-W Align out of range.")
        move = p[at(i+j, offset ÷ 2)]
        if move == 1
            i -= 1; j -= 1; push!(al0, s1[i+1]); push!(al1, s2[j+1])
        elseif move == 2
            j -= 1; push!(al0, 0x00); push!(al1, s2[j+1])
        elseif move == 3
            i -= 1; push!(al0, s1[i+1]); push!(al1, 0x00)
        else
            error("N-W Align out of range.")
        end
    end
    swap && ((al0, al1) = (al1, al0))
    PairAlignment(reverse!(al0), reverse!(al1))
end

"""
Traceback of `nwalign_vectorized` that records only the substitutions of the
query (the second sequence before any swap) against the reference, with the
query positions and transitions of `CompactSub(al2subs-alignment)`.
"""
function compact_traceback(p, s1, s2, len1, len2, ncol, start_col, swap)
    pos = UInt16[]; trans = UInt8[]
    i = len1; j = len2
    @inbounds while i > 0 || j > 0
        offset = 2*start_col + j - i
        offset >= 0 || error("N-W Align out of range.")
        move = p[(i+j)*ncol + offset ÷ 2 + 1]
        if move == 1
            i -= 1; j -= 1
            a = s1[i+1]; b = s2[j+1]
            if a != b  # reference base → query base, at the query position
                ref, query, at_query = swap ? (b, a, i + 1) : (a, b, j + 1)
                push!(pos, at_query); push!(trans, 4*(ref-1) + query)
            end
        elseif move == 2
            j -= 1
        elseif move == 3
            i -= 1
        else
            error("N-W Align out of range.")
        end
    end
    CompactSub(reverse!(pos), reverse!(trans))
end

# ---------------------------------------------------------------------------
# K-mers (kmers.cpp). Counts are dense vectors of 4^5 bytes, as upstream's
# `kmer8`; the rare counts above 255 are kept exactly in a short list, so the
# shared count needs no 16-bit fallback. `kord` keeps the k-mers in order.

function kmer_words(seq::Vector{UInt8}, k::Int=KMER_SIZE)
    n = length(seq) - k + 1
    words = Vector{UInt16}(undef, max(n, 0))
    @inbounds for i in 1:n
        w = 0
        for j in i:i+k-1
            w = 4w + seq[j] - 1
        end
        words[i] = w % UInt16
    end
    words
end

const N_KMERS = 4^KMER_SIZE

"Write k-mer counts saturated at 255 into `dense`; return the exact counts above 255."
function kmer_counts!(dense::AbstractVector{UInt8}, words::Vector{UInt16})
    extra = Dict{UInt16,Int}()
    @inbounds for w in words
        if dense[w+1] < 0xff
            dense[w+1] += 0x01
        else
            extra[w] = get(extra, w, 0) + 1
        end
    end
    sort!([w => 255 + c for (w, c) in extra]; by=first)
end

# Sum over k-mers of min(a, b), in blocks of 256 whose byte minima cannot
# overflow a 16-bit lane (256 × 255 < 2^16): vectorizes as upstream's SSE.
# Fixed trip counts: with a range computed from the length LLVM did not vectorize.
function dense_min_sum(a::AbstractVector{UInt8}, b::AbstractVector{UInt8})
    length(a) == length(b) == N_KMERS || throw(DimensionMismatch("k-mer count vectors"))
    total = 0
    @inbounds for block in 0:N_KMERS÷256-1
        acc = UInt16(0); offset = 256block
        @simd for i in 1:256
            acc += UInt16(min(a[offset+i], b[offset+i]))
        end
        total += Int(acc)
    end
    total
end

function true_count(dense::AbstractVector{UInt8}, overflow, k::Integer)
    c = dense[k+1]
    c < 0xff && return Int(c)
    for (word, count) in overflow
        word == k && return count
    end
    255
end

"Sorted multiset intersection; reference for the dense shared count."
function kmer_dot(a::Vector{UInt16}, b::Vector{UInt16})
    i = j = 1; common = 0
    @inbounds while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            common += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    common
end

"""
Ordered k-mer matches over the shorter sequence (upstream `kord_dist_SSEi`),
from the sequences: the k-mers at a position are equal exactly when the last
`k` bases ending there are all equal.
"""
function kord_dot(a::Vector{UInt8}, b::Vector{UInt8}, k::Int=KMER_SIZE)
    run = 0; common = 0
    @inbounds for i in 1:min(length(a), length(b))
        run = a[i] == b[i] ? run + 1 : 0
        common += run >= k
    end
    common
end

"Ordered k-mer matches of two k-mer vectors (reference for `kord_dot` on sequences)."
function kord_dot(a::Vector{UInt16}, b::Vector{UInt16})
    common = 0
    @inbounds for i in 1:min(length(a), length(b))
        common += a[i] == b[i]
    end
    common
end

# ---------------------------------------------------------------------------
# Substitutions (al2subs)

"Substitutions of a query against a reference; `map` is reference → query position (0 = gap)."
struct Sub
    map::Vector{UInt16}
    pos::Vector{UInt16}  # reference positions
    nt0::Vector{UInt8}
    nt1::Vector{UInt8}
end
nsubs(s::Sub) = length(s.pos)

function al2subs(al::PairAlignment)
    map = UInt16[]; pos = UInt16[]; nt0 = UInt8[]; nt1 = UInt8[]
    i0 = 0; i1 = 0
    @inbounds for k in eachindex(al.al0)
        a = al.al0[k]; b = al.al1[k]
        a != 0 && (i0 += 1)
        b != 0 && (i1 += 1)
        a != 0 && push!(map, b != 0 ? i1 : 0)
        if a != 0 && b != 0 && a != b
            push!(pos, i0); push!(nt0, a); push!(nt1, b)
        end
    end
    Sub(map, pos, nt0, nt1)
end

"Compact substitutions for lambda: query positions and 16-row transition indices."
struct CompactSub
    pos::Vector{UInt16}
    trans::Vector{UInt8}
end
CompactSub(s::Sub) = CompactSub(UInt16[s.map[p] for p in s.pos],
                                UInt8[4(s.nt0[k]-1) + s.nt1[k] for k in eachindex(s.pos)])

function CompactSub(al::PairAlignment)
    n = 0
    @inbounds for k in eachindex(al.al0)
        a = al.al0[k]; b = al.al1[k]
        n += a != 0 && b != 0 && a != b
    end
    pos = Vector{UInt16}(undef, n); trans = Vector{UInt8}(undef, n)
    i1 = 0; s = 0
    @inbounds for k in eachindex(al.al0)
        a = al.al0[k]; b = al.al1[k]
        b != 0 && (i1 += 1)
        if a != 0 && b != 0 && a != b
            s += 1; pos[s] = i1; trans[s] = 4(a-1) + b
        end
    end
    CompactSub(pos, trans)
end

struct AlignOptions
    match::Int
    mismatch::Int
    gap::Int
    homo_gap::Int
    band::Int
    use_kmers::Bool
    kdist_cutoff::Float64
    vectorized::Bool
    gapless::Bool
end

"Column of the shared k-mer count matrix (contiguous, so the scans vectorize)."
const KmerColumn = SubArray{UInt8,1,Matrix{UInt8},Tuple{Base.Slice{Base.OneTo{Int}},Int},true}

"Unique sequence prepared for comparisons."
struct RawSeq
    seq::Vector{UInt8}
    qual::Vector{UInt8}
    kmer::KmerColumn                    # dense counts, saturated at 255
    overflow::Vector{Pair{UInt16,Int}}  # exact counts above 255
end

"""
Uniques prepared for comparisons. Their k-mer counts are the columns of one
4^5 × n matrix: scanning all uniques against a centre then streams through
contiguous memory, about twice as fast as separate vectors.
"""
function raw_seqs(seqs::Vector{Vector{UInt8}}, quals::Vector{Vector{UInt8}}, use_kmers::Bool;
                  threaded::Bool=true)
    n = length(seqs)
    counts = zeros(UInt8, use_kmers ? N_KMERS : 0, n)
    raws = Vector{RawSeq}(undef, n)
    build(i) = begin
        column = view(counts, :, i)
        overflow = use_kmers ? kmer_counts!(column, kmer_words(seqs[i])) : Pair{UInt16,Int}[]
        raws[i] = RawSeq(seqs[i], quals[i], column, overflow)
    end
    if threaded && Threads.nthreads() > 1 && n >= 256
        Threads.@threads for i in 1:n
            build(i)
        end
    else
        foreach(build, 1:n)
    end
    raws
end
RawSeq(seq::Vector{UInt8}, qual::Vector{UInt8}, use_kmers::Bool) = only(raw_seqs([seq], [qual], use_kmers))

"Shared k-mer count of two uniques (sum of minimum multiplicities)."
function kmer_dot(a::RawSeq, b::RawSeq)
    total = dense_min_sum(a.kmer, b.kmer)
    (isempty(a.overflow) && isempty(b.overflow)) && return total
    # Saturated entries: replace the byte minimum by the exact one.
    for (k, _) in a.overflow
        total += min(true_count(a.kmer, a.overflow, k), true_count(b.kmer, b.overflow, k)) -
                 Int(min(a.kmer[k+1], b.kmer[k+1]))
    end
    for (k, _) in b.overflow
        any(p -> p.first == k, a.overflow) && continue
        total += min(true_count(a.kmer, a.overflow, k), true_count(b.kmer, b.overflow, k)) -
                 Int(min(a.kmer[k+1], b.kmer[k+1]))
    end
    total
end

"Upstream `raw_align`: k-mer screen, gapless shortcut, then the chosen aligner."
function raw_align(r0::RawSeq, r1::RawSeq, o::AlignOptions, use_kmers::Bool, kdist_cutoff::Float64,
                   scratch::AlignScratch=AlignScratch(); compact::Bool=false)
    kdist = 0.0; gapless = false
    if use_kmers
        denominator = min(length(r0.seq), length(r1.seq)) - KMER_SIZE + 1.0
        dot = kmer_dot(r0, r1)
        kdist = 1.0 - dot/denominator
        # The ordered count only matters for pairs that pass the screen.
        o.gapless && kdist <= kdist_cutoff && (gapless = kord_dot(r0.seq, r1.seq) == dot)
    end
    if use_kmers && kdist > kdist_cutoff
        nothing
    elseif o.band == 0 || gapless
        nwalign_gapless(r0.seq, r1.seq)
    elseif o.vectorized
        compact ? nwalign_vectorized(r0.seq, r1.seq, o.match, o.mismatch, o.gap, 0, o.band; scratch,
                                     compact=Val(true)) :
                  nwalign_vectorized(r0.seq, r1.seq, o.match, o.mismatch, o.gap, 0, o.band; scratch)
    else
        nwalign_endsfree(r0.seq, r1.seq, o.match, o.mismatch, o.gap, o.band;
                         homo_gap=o.homo_gap != o.gap && o.homo_gap <= 0 ? o.homo_gap : o.gap)
    end
end

function sub_new(r0::RawSeq, r1::RawSeq, o::AlignOptions, use_kmers::Bool, kdist_cutoff::Float64,
                 scratch::AlignScratch=AlignScratch())
    al = raw_align(r0, r1, o, use_kmers, kdist_cutoff, scratch)
    al === nothing ? nothing : al2subs(al)
end

function compact_sub_new(r0::RawSeq, r1::RawSeq, o::AlignOptions, use_kmers::Bool,
                         kdist_cutoff::Float64, scratch::AlignScratch)
    al = raw_align(r0, r1, o, use_kmers, kdist_cutoff, scratch; compact=true)
    al === nothing ? nothing : al isa CompactSub ? al : CompactSub(al)
end

"Index ranges processed by one task each, in increasing order."
work_chunks(n::Int, threaded::Bool) = threaded && Threads.nthreads() > 1 && n >= 64 ?
    collect(Iterators.partition(1:n, cld(n, 4*Threads.nthreads()))) : [1:n]

"Run `f(scratch, i)` for `i in 1:n`, in chunks with one scratch per chunk."
function foreach_scratch(f, n::Int, threaded::Bool)
    if threaded && Threads.nthreads() > 1 && n >= 64
        chunks = collect(Iterators.partition(1:n, cld(n, 4*Threads.nthreads())))
        Threads.@threads for chunk in chunks
            local scratch = AlignScratch()
            for i in chunk
                f(scratch, i)
            end
        end
    else
        local serial = AlignScratch()
        for i in 1:n
            f(serial, i)
        end
    end
end

# ---------------------------------------------------------------------------
# Lambda and abundance p-values (pval.cpp)

"Probability that the center produces the raw: product over raw positions, in order."
function compute_lambda(raw::RawSeq, sub::CompactSub, err::Matrix{Float64})
    lambda = 1.0; s = 1; ns = length(sub.pos)
    @inbounds for p in eachindex(raw.seq)
        if s <= ns && sub.pos[s] == p
            t = Int(sub.trans[s]); s += 1
        else
            t = 5*(raw.seq[p]-1) + 1
        end
        lambda *= err[t, raw.qual[p]+1]
    end
    0 <= lambda <= 1 || error("Bad lambda.")
    lambda
end

"Upstream `calc_pA`: P(X ≥ reads), conditional on presence unless `prior`."
function calc_pA(reads::Integer, expected::Float64, prior::Bool)
    pval = ccdf(Poisson(expected), reads - 1)
    # As R's pgamma, results close to underflow are recomputed in log space,
    # keeping subnormal values that order the most significant sequences.
    pval < floatmin(Float64)/eps(Float64) && (pval = exp(logccdf(Poisson(expected), reads - 1)))
    if !prior
        norm = 1.0 - exp(-expected)
        norm < TAIL_APPROX_CUTOFF && (norm = expected - 0.5*expected*expected)
        pval = pval/norm
    end
    pval
end

# ---------------------------------------------------------------------------
# Partition (containers.cpp, cluster.cpp)

struct Comparison
    i::Int32        # cluster
    index::Int32    # raw
    lambda::Float64
    hamming::Int32  # typemax for a screened (NULL) comparison
end
const NULL_HAMMING = typemax(Int32)
const NO_COMPARISON = Comparison(0, 0, 0.0, 0)

mutable struct Cluster
    center::Int
    raws::Vector{Int}
    reads::Int
    update_e::Bool
    check_locks::Bool
    birth_from::Int
    birth_pval::Float64
    birth_fold::Float64
    birth_e::Float64
    birth_comp::Comparison
    comps::Vector{Comparison}
end
Cluster() = Cluster(0, Int[], 0, true, true, 0, 0.0, 1.0, 0.0, NO_COMPARISON, Comparison[])

"""
Substitutions of the uniques compared with one centre (`nothing`: screened
out by k-mers), with indices in increasing order.
"""
struct CenterComparisons
    index::Vector{Int32}
    subs::Vector{Union{Nothing,CompactSub}}
end
Base.length(c::CenterComparisons) = length(c.index)

"Merge two index-ordered comparison lists (indices are disjoint)."
function merge_comparisons(a::CenterComparisons, index::Vector{Int32}, subs)
    isempty(index) && return a
    n = length(a.index) + length(index)
    out = CenterComparisons(Vector{Int32}(undef, n), Vector{Union{Nothing,CompactSub}}(undef, n))
    i = j = 1
    for k in 1:n
        if j > length(index) || (i <= length(a.index) && a.index[i] < index[j])
            out.index[k] = a.index[i]; out.subs[k] = a.subs[i]; i += 1
        else
            out.index[k] = index[j]; out.subs[k] = subs[j]; j += 1
        end
    end
    out
end

"Mean qualities of unique `i` rounded half away from zero, as upstream `raw_new`."
rounded_qualities(q::MeanQualities, i::Int) = q.counts[i] == 1 ? UInt8.(q.sums[i]) :
    UInt8[round(Int, s / q.counts[i], RoundNearestTiesAway) for s in q.sums[i]]
rounded_qualities(q::AbstractVector, i::Int) = UInt8[round(Int, x, RoundNearestTiesAway) for x in q[i]]

"Highest mean quality of a dereplicated sample."
max_quality(q::MeanQualities) =
    maximum((maximum(q.sums[i]; init=0) / q.counts[i] for i in eachindex(q.sums)); init=0.0)
max_quality(q::AbstractVector) = maximum(x -> maximum(x; init=0.0), q; init=0.0)

"Caller-owned sequence data and alignment caches, reusable across error-learning rounds."
mutable struct InferenceCache
    raws::Vector{RawSeq}
    order::Vector{Int}                       # processing index → Derep index
    options::Union{Nothing,AlignOptions}
    cache_alignments::Bool
    compare::Vector{Union{Nothing,CenterComparisons}}  # by centre
    final::Vector{Union{Nothing,Tuple{Int32,Sub}}}     # by unique: (centre, alignment)
end

function InferenceCache(d::Derep; useKmers=true, cacheAlignments=true, threaded=true, kwargs...)
    n = length(d.sequences)
    # Process uniques by decreasing abundance, ties in input order (derepFastq).
    order = sortperm(1:n; by=i -> (-d.counts[i], i))
    raws = raw_seqs([encode_dna(d.sequences[i]) for i in order],
                    [rounded_qualities(d.qualities, i) for i in order],
                    useKmers; threaded)
    InferenceCache(raws, order, nothing, cacheAlignments,
                   Vector{Union{Nothing,CenterComparisons}}(nothing, n),
                   Vector{Union{Nothing,Tuple{Int32,Sub}}}(nothing, n))
end

function check_cache!(cache::InferenceCache, d::Derep, o::AlignOptions, cache_alignments::Bool)
    length(cache.raws) == length(d.sequences) && cache.cache_alignments == cache_alignments &&
        (!o.use_kmers || all(r -> length(r.kmer) == N_KMERS, cache.raws)) ||
        throw(ArgumentError("Invalid inference cache"))
    if cache.options === nothing
        cache.options = o
    elseif cache.options != o
        # Alignments depend on the options: start again rather than mix them.
        fill!(cache.compare, nothing); fill!(cache.final, nothing); cache.options = o
    end
    cache
end

mutable struct Partition
    raws::Vector{RawSeq}
    reads::Vector{Int}
    prior::Vector{Bool}
    p::Vector{Float64}
    E_minmax::Vector{Float64}
    comp::Vector{Comparison}
    lock::Vector{Bool}
    correct::Vector{Bool}
    clusters::Vector{Cluster}
    total::Int
    nalign::Int
    omegaA::Float64
    omegaP::Float64
end

function add_raw!(c::Cluster, b::Partition, r::Int)
    push!(c.raws, r); c.reads += b.reads[r]; c.update_e = true
    c
end

function pop_raw!(c::Cluster, b::Partition, k::Int)
    r = c.raws[k]
    c.raws[k] = c.raws[end]; pop!(c.raws)
    c.reads -= b.reads[r]; c.update_e = true
    r
end

function assign_center!(c::Cluster, b::Partition)
    c.center = 0; max_reads = 0
    for r in c.raws
        b.lock[r] = false
        if b.reads[r] > max_reads
            c.center = r; max_reads = b.reads[r]
        end
    end
    c.check_locks = true
    c
end

function Partition(raws, reads, prior, omegaA, omegaP)
    n = length(raws)
    b = Partition(raws, reads, prior, zeros(n), fill(-999.0, n), fill(NO_COMPARISON, n),
                  falses(n), trues(n), Cluster[], sum(reads), 0, omegaA, omegaP)
    c = Cluster()
    c.birth_e = b.total
    for r in 1:n
        add_raw!(c, b, r)
    end
    push!(b.clusters, c)
    assign_center!(c, b)
    b
end

struct DadaOptions
    align::AlignOptions
    err::Matrix{Float64}
    greedy::Bool
    detect_singletons::Bool
    min_fold::Float64
    min_hamming::Int
    min_abund::Int
    threaded::Bool
end

function compare_one(b::Partition, o::DadaOptions, center::Int, index::Int,
                     cutoff::Float64, scratch::AlignScratch)
    center_reads = b.reads[center]
    if o.greedy && (b.reads[index] > center_reads || b.lock[index])
        return 0.0, NULL_HAMMING, false, nothing
    end
    sub = compact_sub_new(b.raws[center], b.raws[index], o.align, o.align.use_kmers, cutoff, scratch)
    sub === nothing && return 0.0, NULL_HAMMING, true, sub
    compute_lambda(b.raws[index], sub, o.err), Int32(length(sub.pos)), true, sub
end

"Per-task output of `compare!`: stored comparisons and newly computed substitutions."
struct CompareChunk
    comps::Vector{Comparison}
    fresh_index::Vector{Int32}
    fresh_subs::Vector{Union{Nothing,CompactSub}}
end
CompareChunk() = CompareChunk(Comparison[], Int32[], Union{Nothing,CompactSub}[])

"""
Upstream `b_compare_parallel`: compare every raw with the center of cluster `i`.
Each task also applies upstream's storing rule to its own indices
(`E_minmax` and `comp` of an index are only touched by that index), and the
stored comparisons are concatenated in chunk order, that is by index as in R.
"""
function compare!(b::Partition, i::Int, cache::InferenceCache, o::DadaOptions, cutoff::Float64)
    n = length(b.raws); c = b.clusters[i]; center = c.center
    cached = cache.cache_alignments ? cache.compare[center] : nothing
    chunks = work_chunks(n, o.threaded)
    out = [CompareChunk() for _ in chunks]
    work(k) = begin
        local scratch = AlignScratch()
        chunk = chunks[k]; result = out[k]
        # Cursor into the cached, index-ordered substitutions of this centre.
        cursor = cached === nothing ? 1 : searchsortedfirst(cached.index, Int32(first(chunk)))
        center_reads = b.reads[center]
        for index in chunk
            known = nothing; hit = false
            if cached !== nothing
                while cursor <= length(cached.index) && cached.index[cursor] < index
                    cursor += 1
                end
                if cursor <= length(cached.index) && cached.index[cursor] == index
                    hit = true; known = cached.subs[cursor]
                end
            end
            lambda, hamming, fresh, sub = hit ?
                compare_cached(b, o, index, center_reads, known) :
                compare_one(b, o, center, index, cutoff, scratch)
            0 <= lambda <= 1 || error("Lambda out-of-range error.")
            if fresh && cache.cache_alignments
                push!(result.fresh_index, Int32(index)); push!(result.fresh_subs, sub)
            end
            if lambda * b.total > b.E_minmax[index]
                e = lambda * center_reads
                e > b.E_minmax[index] && (b.E_minmax[index] = e)
                comp = Comparison(i, index, lambda, hamming)
                push!(result.comps, comp)
                (i == 1 || index == center) && (b.comp[index] = comp)
            end
        end
    end
    if length(chunks) > 1
        Threads.@threads for k in eachindex(chunks)
            work(k)
        end
    else
        work(1)
    end
    b.nalign += n
    for result in out
        append!(c.comps, result.comps)
    end
    if cache.cache_alignments
        index = reduce(vcat, (r.fresh_index for r in out); init=Int32[])
        subs = reduce(vcat, (r.fresh_subs for r in out); init=Union{Nothing,CompactSub}[])
        cache.compare[center] = cached === nothing ? CenterComparisons(index, subs) :
                                merge_comparisons(cached, index, subs)
    end
    b
end

"Comparison from a cached substitution (greedy skips are checked first, as upstream)."
function compare_cached(b::Partition, o::DadaOptions, index::Int, center_reads::Int, sub)
    if o.greedy && (b.reads[index] > center_reads || b.lock[index])
        return 0.0, NULL_HAMMING, false, nothing
    end
    sub === nothing && return 0.0, NULL_HAMMING, false, nothing
    compute_lambda(b.raws[index], sub, o.err), Int32(length(sub.pos)), false, sub
end

"Upstream `b_shuffle2`: move raws to the cluster expected to produce the most reads."
function shuffle!(b::Partition)
    n = length(b.raws)
    first = b.clusters[1].comps
    compmax = Vector{Comparison}(undef, n); emax = zeros(n)
    for index in 1:n
        comp = first[index]
        comp.index == index || error("Missing initial comparison")
        compmax[index] = comp; emax[index] = comp.lambda * b.clusters[1].reads
    end
    for i in 2:length(b.clusters)
        reads = b.clusters[i].reads
        for comp in b.clusters[i].comps
            e = comp.lambda * reads
            if e > emax[comp.index]
                compmax[comp.index] = comp; emax[comp.index] = e
            end
        end
    end
    shuffled = false
    for i in eachindex(b.clusters)
        c = b.clusters[i]
        for k in length(c.raws):-1:1
            r = c.raws[k]
            target = Int(compmax[r].i)
            target == i && continue
            r == c.center && continue
            pop_raw!(c, b, k)
            add_raw!(b.clusters[target], b, r)
            b.comp[r] = compmax[r]
            shuffled = true
        end
    end
    shuffled
end

function get_pA(b::Partition, r::Int, c::Cluster, detect_singletons::Bool)
    comp = b.comp[r]
    if b.reads[r] == 1 && !b.prior[r] && !detect_singletons
        1.0
    elseif comp.hamming == 0
        1.0
    elseif comp.lambda == 0
        0.0
    else
        calc_pA(b.reads[r], comp.lambda * c.reads, b.prior[r] || detect_singletons)
    end
end

"Upstream `b_p_update`, including greedy locking on the first update of a cluster."
function p_update!(b::Partition, o::DadaOptions)
    for c in b.clusters
        if c.update_e
            for r in c.raws
                b.p[r] = get_pA(b, r, c, o.detect_singletons)
            end
            c.update_e = false
        end
        if o.greedy && c.check_locks
            center_reads = b.reads[c.center]
            for r in c.raws
                center_reads * b.comp[r].lambda > b.reads[r] && (b.lock[r] = true)
                r == c.center && (b.lock[r] = true)
            end
            c.check_locks = false
        end
    end
    b
end

"Most significant raws (abundance and prior), as searched by upstream `b_bud`."
function bud_candidates(b::Partition, o::DadaOptions)
    minraw = minprior = b.clusters[1].center
    mini = minr = mini_prior = minr_prior = 0
    for (i, c) in enumerate(b.clusters), k in 2:length(c.raws)  # position 1 is the center
        r = c.raws[k]
        b.reads[r] < o.min_abund && continue
        comp = b.comp[r]
        comp.hamming >= o.min_hamming || continue
        o.min_fold <= 1 || b.reads[r] >= o.min_fold * comp.lambda * c.reads || continue
        if b.p[r] < b.p[minraw] || (b.p[r] == b.p[minraw] && b.reads[r] > b.reads[minraw])
            mini = i; minr = k; minraw = r
        end
        if b.prior[r] && (b.p[r] < b.p[minprior] ||
                          (b.p[r] == b.p[minprior] && b.reads[r] > b.reads[minprior]))
            mini_prior = i; minr_prior = k; minprior = r
        end
    end
    pA = b.p[minraw] * length(b.raws)
    pP = b.p[minprior]
    if pA < b.omegaA && mini > 0
        (mini, minr, minraw, pA)
    elseif pP < b.omegaP && mini_prior > 0
        (mini_prior, minr_prior, minprior, pP)
    else
        nothing
    end
end

"Upstream `b_bud`: split the most significant raw into a new cluster; 0 if none."
function bud!(b::Partition, o::DadaOptions)
    found = bud_candidates(b, o)
    found === nothing && return 0
    i, k, r, pval = found
    comp = b.comp[r]
    expected = comp.lambda * b.clusters[i].reads
    pop_raw!(b.clusters[i], b, k)
    c = Cluster()
    c.birth_from = i; c.birth_pval = pval; c.birth_fold = b.reads[r]/expected
    c.birth_e = expected; c.birth_comp = comp
    push!(b.clusters, c)
    add_raw!(c, b, r)
    assign_center!(c, b)
    length(b.clusters)
end

"Upstream `run_dada`; returns the partition and whether it converged."
function run_dada!(b::Partition, cache::InferenceCache, o::DadaOptions, max_clust::Int)
    # Everyone is aligned with the initial center, without the k-mer screen.
    compare!(b, 1, cache, o, 1.0)
    p_update!(b, o)
    max_clust < 1 && (max_clust = length(b.raws))
    converged = true
    while length(b.clusters) < max_clust
        newi = bud!(b, o)
        newi == 0 && break
        compare!(b, newi, cache, o, o.align.kdist_cutoff)
        nshuffle = 0
        while true
            shuffled = shuffle!(b)
            shuffled && (nshuffle += 1) < MAX_SHUFFLE && continue
            shuffled && (converged = false)  # Reached the maximum number of shuffles.
            break
        end
        p_update!(b, o)
    end
    length(b.clusters) >= max_clust && bud_candidates(b, o) !== nothing && (converged = false)
    converged
end

"Alignment of a raw with its final center, without the k-mer screen (cached)."
function final_sub(b::Partition, cache::InferenceCache, o::DadaOptions, center::Int, r::Int,
                   scratch::AlignScratch=AlignScratch())
    # Each unique has one final centre per run: its slot is written by one task.
    if cache.cache_alignments
        entry = cache.final[r]
        entry !== nothing && entry[1] == center && return entry[2]
    end
    sub = sub_new(b.raws[center], b.raws[r], o.align, false, 1.0, scratch)
    cache.cache_alignments && (cache.final[r] = (Int32(center), sub))
    sub
end
