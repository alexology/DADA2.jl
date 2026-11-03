const DNA = UInt8[0x41, 0x43, 0x47, 0x54]
@inline function baseindex(b::UInt8)
    b == 0x41 && return 1
    b == 0x43 && return 2
    b == 0x47 && return 3
    b == 0x54 && return 4
    throw(ArgumentError("Expected A/C/G/T, found $(Char(b))"))
end
const IUPAC_DNA = codeunits("ACGTRYSWKMBDHVN")
const DNAInput = Union{AbstractString,BioSequences.NucSeq}
@inline ascii_upper(b::UInt8) = 0x61 <= b <= 0x7a ? b - 0x20 : b
# Class of each raw byte, after upper-casing: 0x01 = A/C/G/T, 0x02 = N,
# 0x04 = IUPAC code; 0x80 marks a lower-case letter. Non-ASCII bytes have no
# class, so ASCII upper-casing is enough.
const DNA_CLASS = let t = zeros(UInt8, 256)
    for b in 0x00:0xff
        u = ascii_upper(b)
        u in DNA && (t[b+1] |= 0x01)
        u == 0x4e && (t[b+1] |= 0x02)
        u in IUPAC_DNA && (t[b+1] |= 0x04)
        u != b && (t[b+1] |= 0x80)
    end
    t
end
dna_allowed(ambiguous, iupac) = 0x01 | (ambiguous ? 0x02 : 0x00) | (iupac ? 0x04 : 0x00)
# Upper-cased sequence, or `nothing` if a byte is not allowed.
function checked_dna(seq::String, allowed::UInt8)
    cu = codeunits(seq)
    seen = 0x00; bad = false
    @inbounds @simd for i in eachindex(cu)
        c = DNA_CLASS[cu[i]+1]
        seen |= c
        bad |= c & allowed == 0x00
    end
    bad ? nothing : seen & 0x80 != 0x00 ? String(map(ascii_upper, cu)) : seq
end
function validated_dna(s::AbstractString; ambiguous=false, iupac=false)
    seq = checked_dna(String(s), dna_allowed(ambiguous, iupac))
    seq === nothing && throw(ArgumentError("Sequence must contain A/C/G/T$(ambiguous ? "/N" : "")"))
    seq
end
validated_dna(s::BioSequences.NucSeq;kwargs...)=validated_dna(String(s);kwargs...)

struct FastqRecord
    id::String
    sequence::String
    quality::Vector{UInt8}
    function FastqRecord(id::AbstractString, seq::AbstractString, quality::AbstractVector{<:Integer})
        s = validated_dna(seq; ambiguous=true)
        ncodeunits(s) == length(quality) || throw(ArgumentError("Sequence/quality length mismatch"))
        all(q -> 0 <= q <= 93, quality) || throw(ArgumentError("Phred scores must be in 0:93"))
        new(String(id), s, UInt8.(quality))
    end
end
FastqRecord(id::AbstractString,seq::BioSequences.NucSeq,quality::AbstractVector{<:Integer})=
    FastqRecord(id,String(seq),quality)

"Dereplicated reads. `map` maps original read order to unique sequence indices."
struct Derep
    sequences::Vector{String}
    counts::Vector{Int}
    qualities::Vector{Vector{Float64}}
    map::Vector{Int}
    function Derep(sequences, counts, qualities, map=Int[])
        length(sequences) == length(counts) == length(qualities) || throw(ArgumentError("Inconsistent derep dimensions"))
        seq = [validated_dna(s) for s in sequences]
        length(unique(seq)) == length(seq) || throw(ArgumentError("Duplicate unique sequences"))
        all(>(0), counts) || throw(ArgumentError("Counts must be positive"))
        qs = [Float64.(q) for q in qualities]
        all(i -> length(qs[i]) == ncodeunits(seq[i]), eachindex(seq)) || throw(ArgumentError("Quality lengths differ"))
        all(q -> all(x -> isfinite(x) && 0 <= x <= 93, q), qs) || throw(ArgumentError("Invalid qualities"))
        all(i -> 1 <= i <= length(seq), map) || throw(ArgumentError("Invalid read map"))
        if !isempty(map)
            observed = zeros(Int, length(seq))
            for i in map
                observed[i] += 1
            end
            observed == counts || throw(ArgumentError("Map disagrees with counts"))
        end
        new(seq, Int.(counts), qs, Int.(map))
    end
end

"16 × (maxQ+1) transition probabilities: A→A, A→C, …, T→T."
struct ErrorModel
    rates::Matrix{Float64}
    logrates::Matrix{Float64}
    function ErrorModel(rates::AbstractMatrix{<:Real})
        size(rates, 1) == 16 && size(rates, 2) > 0 || throw(ArgumentError("Expected 16 × Q matrix"))
        r = Float64.(rates)
        all(x -> isfinite(x) && 0 <= x <= 1, r) || throw(ArgumentError("Invalid error probabilities"))
        for q in axes(r, 2), b in 0:3
            isapprox(sum(@view r[4b+1:4b+4,q]), 1; atol=1e-8) || throw(ArgumentError("Transition rows must sum to one per source base and quality"))
        end
        new(r, log.(r))
    end
    # Upstream intermediate models (all ones, or self-transitions forced to
    # one) do not sum to one; they are only used inside error learning.
    ErrorModel(rates::Matrix{Float64}, ::Val{:unchecked}) = new(rates, log.(rates))
end

"""
Diagnostics of each ASV, as R's `dada(...)\$clustering`, in the order of
`DadaResult.sequences`: `abundance`, reads of uniques without (`n0`) or with
one (`n1`) substitution from the centre, `nunq`, post-hoc `pval`, and the
birth of the partition (`birth_from` index in this result, `birth_pval`,
`birth_fold`, `birth_ham`, `birth_qave`), missing for the first partition.
"""
const ClusteringTable = NamedTuple{(:sequence, :abundance, :n0, :n1, :nunq, :pval, :birth_from,
    :birth_pval, :birth_fold, :birth_ham, :birth_qave),
    Tuple{Vector{String}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Float64},
          Vector{Union{Missing,Int}}, Vector{Union{Missing,Float64}}, Vector{Union{Missing,Float64}},
          Vector{Union{Missing,Int}}, Vector{Union{Missing,Float64}}}}

empty_clustering() = ClusteringTable((String[], Int[], Int[], Int[], Int[], Float64[],
    Union{Missing,Int}[], Union{Missing,Float64}[], Union{Missing,Float64}[],
    Union{Missing,Int}[], Union{Missing,Float64}[]))

"Rows `rows` of a clustering table, with `birth_from` renumbered by `index` (0 = absent)."
function select_clustering(t::ClusteringTable, rows::AbstractVector{Int}, index::AbstractVector{Int})
    birth = Union{Missing,Int}[x === missing || index[x] == 0 ? missing : index[x] for x in t.birth_from[rows]]
    ClusteringTable((t.sequence[rows], t.abundance[rows], t.n0[rows], t.n1[rows], t.nunq[rows],
        t.pval[rows], birth, t.birth_pval[rows], t.birth_fold[rows], t.birth_ham[rows],
        t.birth_qave[rows]))
end

struct DadaResult
    sequences::Vector{String}
    counts::Vector{Int}
    assignments::Vector{Int}
    derep::Derep
    transitions::Matrix{Int}
    errors::ErrorModel
    converged::Bool
    comparisons::Int
    rejected::Int
    clustering::ClusteringTable
end

struct SequenceTable
    samples::Vector{String}
    sequences::Vector{String}
    counts::Matrix{Int}
    function SequenceTable(samples, sequences, counts)
        size(counts) == (length(samples), length(sequences)) || throw(ArgumentError("Table dimensions disagree"))
        length(unique(samples)) == length(samples) || throw(ArgumentError("Duplicate sample names"))
        length(unique(sequences)) == length(sequences) || throw(ArgumentError("Duplicate table sequences"))
        all(>=(0), counts) || throw(ArgumentError("Negative abundance"))
        new(String.(samples), String.(sequences), Int.(counts))
    end
end

"Pairs of forward/reverse ASVs, as the data frame of R `mergePairs(..., returnRejects=TRUE)`."
const MergeTable = NamedTuple{(:sequence, :abundance, :forward, :reverse, :nmatch, :nmismatch,
                               :nindel, :prefer, :accept),
    Tuple{Vector{String}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int},
          Vector{Int}, Vector{Union{Missing,Int}}, Vector{Bool}}}

"""
Merged reads: `sequences`/`counts` of accepted merges (identical sequences
summed), read totals `accepted`/`rejected` (rejected includes reads whose
forward or reverse read was not assigned) and every ASV pair in `pairs`,
sorted by decreasing abundance (unmerged pairs have an empty sequence).
"""
struct MergeResult
    sequences::Vector{String}
    counts::Vector{Int}
    accepted::Int
    rejected::Int
    pairs::MergeTable
end
MergeResult(sequences, counts, accepted, rejected) =
    MergeResult(sequences, counts, accepted, rejected,
                MergeTable((String[], Int[], Int[], Int[], Int[], Int[], Int[], Union{Missing,Int}[], Bool[])))

Base.@kwdef struct FilterOptions
    trimLeft::Int = 0
    truncLen::Int = 0
    truncQ::Int = 2
    maxN::Int = 0
    maxEE::Float64 = 2.0
    minLen::Int = 20
    maxLen::Int = typemax(Int)
    rmPhix::Bool = true  # discard reads matching phiX (R `rm.phix`)
end
struct FilterStats
    reads_in::Int
    reads_out::Int
end

struct TaxonomyResult
    sequences::Vector{String}
    taxa::Matrix{Union{Missing,String}}
    bootstrap::Matrix{Int}
end
