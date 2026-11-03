function makeSequenceTable(results::AbstractVector; samples=["sample$i" for i in eachindex(results)])
    length(samples) == length(results) || throw(ArgumentError("Wrong number of sample names"))
    seq = sort!(unique([s for r in results for s in r.sequences]))
    index = Dict(s => i for (i,s) in enumerate(seq))
    counts = zeros(Int,length(results),length(seq))
    for (i,r) in enumerate(results), (s,n) in zip(r.sequences,r.counts)
        counts[i,index[s]] += n
    end
    order = sortperm(eachindex(seq); by=j -> (-sum(@view(counts[:,j])),seq[j]))
    SequenceTable(samples,seq[order],counts[:,order])
end

isspace_byte(b::UInt8) = b == 0x20 || 0x09 <= b <= 0x0d
# Whitespace inside sequence lines is removed only when validation fails, so
# that clean sequences are scanned once.
function validated_fasta_dna(s::String)
    seq = checked_dna(s, dna_allowed(false, true))
    seq !== nothing && return seq
    validated_dna(String(filter(!isspace_byte, codeunits(s))); iupac=true)
end
function readFasta(path::AbstractString; bio=false)
    reader = FASTX.FASTAReader(open_sequence(path,"r"); copy=false); ids = String[]; sequences = String[]
    try
        for rec in reader
            id = strip(FASTX.description(rec))
            isempty(id) && throw(ArgumentError("Empty FASTA header"))
            s = validated_fasta_dna(FASTX.sequence(String,rec))
            isempty(s) && throw(ArgumentError("Empty FASTA sequence"))
            push!(ids,id); push!(sequences,s)
        end
    catch e
        e isa ErrorException || e isa EOFError || rethrow()
        throw(ArgumentError("Invalid FASTA $path: $(sprint(showerror, e))"))
    finally
        close(reader)
    end
    (ids=ids,sequences=bio ? BioSequences.LongDNA{4}.(sequences) : sequences)
end
function writeFasta(path::AbstractString, sequences::AbstractVector{<:DNAInput}; ids=["ASV$i" for i in eachindex(sequences)])
    length(ids) == length(sequences) || throw(ArgumentError("Header count differs"))
    isfile(path) && throw(ArgumentError("Output already exists"))
    seq = [validated_dna(s; iupac=true) for s in sequences]
    any(isempty,seq) && throw(ArgumentError("Empty FASTA sequence"))
    # Headers must read back: no line breaks, not empty or blank (DD7).
    all(id -> !occursin('\n',id) && !occursin('\r',id) && !isempty(strip(id)),ids) ||
        throw(ArgumentError("Invalid FASTA header"))
    io = open_sequence(path,"w")
    try
        for (id,s) in zip(ids,seq)
            println(io, '>', id); println(io,s)
        end
    finally
        close(io)
    end
    path
end
