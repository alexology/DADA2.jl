# Port of R dada2 `isPhiX` and `C_matchRef` (filter.R, filter.cpp). The phiX
# genome is the file distributed with dada2 (data/phix_genome.fa).

const PHIX_PATH = joinpath(@__DIR__, "..", "data", "phix_genome.fa")
include_dependency(PHIX_PATH)

function read_phix(path=PHIX_PATH)
    reader = FASTX.FASTAReader(open(path))
    try
        uppercase(FASTX.sequence(String, only(collect(reader))))
    finally
        close(reader)
    end
end
const PHIX_GENOME = read_phix()

# Packed word of `seq[i:i+k-1]` (two bits per base), or `nothing` with a
# non-ACGT base: such words cannot match the reference.
function packed_word(seq::AbstractVector{UInt8}, i::Int, k::Int)
    w = UInt64(0)
    @inbounds for j in i:i+k-1
        b = seq[j]
        c = b == 0x41 ? 0 : b == 0x43 ? 1 : b == 0x47 ? 2 : b == 0x54 ? 3 : -1
        c < 0 && return nothing
        w = (w << 2) | UInt64(c)
    end
    w
end

"Words of a circular reference, as `C_matchRef` (the first `k` bases are appended)."
function reference_words(ref::String, k::Int)
    seq = codeunits(ref * ref[1:k])
    words = Set{UInt64}()
    for i in 1:ncodeunits(ref)
        w = packed_word(seq, i, k)
        w === nothing || push!(words, w)
    end
    words
end

const PHIX_WORDS = Dict{Int,Tuple{Set{UInt64},Set{UInt64}}}()
const PHIX_LOCK = ReentrantLock()
phix_words(k::Int) = lock(PHIX_LOCK) do
    get!(PHIX_WORDS, k) do
        (reference_words(PHIX_GENOME, k), reference_words(reverseComplement(PHIX_GENOME), k))
    end
end

@inline base_code(b::UInt8) = b == 0x41 ? 0 : b == 0x43 ? 1 : b == 0x47 ? 2 : b == 0x54 ? 3 : -1

"""
Upstream `C_matchRef`: reference words found in `seq`. With `non_overlapping`,
a hit moves the window start `j` to `j + k + 1` (upstream adds `k`, then the
loop increment), so the rolling word restarts there.
"""
function match_ref(seq::AbstractVector{UInt8}, words::Set{UInt64}, k::Int, non_overlapping::Bool)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - 1
    hits = 0; w = UInt64(0); run = 0; i = 1; n = length(seq)
    @inbounds while i <= n
        c = base_code(seq[i])
        if c < 0
            run = 0; w = UInt64(0)
        else
            w = ((w << 2) | UInt64(c)) & mask; run += 1
        end
        if run >= k && w in words
            hits += 1
            if non_overlapping
                i += 1  # the next window starts two bases after this one ends
                run = 0; w = UInt64(0)
            end
        end
        i += 1
    end
    hits
end

"""
    isPhiX(sequences; wordSize=16, minMatches=2, nonOverlapping=true)

Port of R `isPhiX`: true for sequences sharing at least `minMatches` words of
`wordSize` bases with the phiX genome, on either strand (counted separately;
with `nonOverlapping`, the search skips past each matching word).
"""
function isPhiX(sequences::AbstractVector{<:DNAInput}; wordSize::Integer=16, minMatches::Integer=2,
                nonOverlapping::Bool=true)
    1 <= wordSize <= 32 || throw(ArgumentError("wordSize must be in 1:32"))
    forward, reverse = phix_words(Int(wordSize))
    map(sequences) do s
        seq = codeunits(uppercase(String(s)))
        match_ref(seq, forward, Int(wordSize), nonOverlapping) >= minMatches ||
            match_ref(seq, reverse, Int(wordSize), nonOverlapping) >= minMatches
    end
end
isPhiX(sequence::DNAInput; kwargs...) = only(isPhiX([sequence]; kwargs...))

# Default screen used by filterAndTrim.
is_phix(seq::AbstractString) = isPhiX(seq)
