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
function validated_dna(s::AbstractString; ambiguous=false, iupac=false)
    seq = uppercase(String(s))
    all(b -> b in DNA || (ambiguous && b == 0x4e) || (iupac && b in IUPAC_DNA), codeunits(seq)) ||
        throw(ArgumentError("Sequence must contain A/C/G/T$(ambiguous ? "/N" : "")"))
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

struct MergeResult
    sequences::Vector{String}
    counts::Vector{Int}
    accepted::Int
    rejected::Int
end

Base.@kwdef struct FilterOptions
    trimLeft::Int = 0
    truncLen::Int = 0
    truncQ::Int = 2
    maxN::Int = 0
    maxEE::Float64 = 2.0
    minLen::Int = 20
    maxLen::Int = typemax(Int)
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
