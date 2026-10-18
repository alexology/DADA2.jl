"Streaming four-line FASTQ reader; caller owns and closes the reader."
mutable struct FastqReader
    io::IO
    line::Int
end
function open_sequence(path::AbstractString, mode::AbstractString)
    raw = open(path, mode)
    if endswith(lowercase(path), ".gz")
        try
            return mode == "r" ? GzipDecompressorStream(raw) : GzipCompressorStream(raw)
        catch
            close(raw)
            rethrow()
        end
    end
    raw
end
FastqReader(path::AbstractString) = FastqReader(open_sequence(path, "r"), 0)
Base.close(r::FastqReader) = close(r.io)
Base.IteratorSize(::Type{FastqReader}) = Base.SizeUnknown()
Base.eltype(::Type{FastqReader}) = FastqRecord
function Base.iterate(r::FastqReader, state=nothing)
    eof(r.io) && return nothing
    lines = String[]
    for _ in 1:4
        eof(r.io) && throw(ArgumentError("Truncated FASTQ at line $(r.line+1)"))
        push!(lines, readline(r.io))
        r.line += 1
    end
    startswith(lines[1], "@") && startswith(lines[3], "+") || throw(ArgumentError("Invalid FASTQ at line $(r.line-3)"))
    qual = [Int(c)-33 for c in codeunits(lines[4])]
    (FastqRecord(lines[1][2:end], lines[2], qual), nothing)
end
function write_record(io::IO, r::FastqRecord)
    println(io, '@', r.id)
    println(io, r.sequence)
    println(io, '+')
    for q in r.quality
        write(io, q + UInt8(33))
    end
    write(io, UInt8(10))
end
const PHRED_PROB = [10.0^(-q/10) for q in 0:93]
expectedErrors(q::AbstractVector{<:Integer}) = sum(x -> PHRED_PROB[Int(x)+1], q; init=0.0)
expectedErrors(r::FastqRecord) = expectedErrors(r.quality)
function check_options(o::FilterOptions)
    o.trimLeft >= 0 && o.truncLen >= 0 && 0 <= o.truncQ <= 93 && o.maxN >= 0 &&
        o.maxEE >= 0 && o.minLen >= 1 && o.maxLen >= o.minLen || throw(ArgumentError("Invalid filter options"))
end
function filtered(r::FastqRecord, o::FilterOptions)
    last = length(r.quality)
    o.truncLen > 0 && last < o.truncLen && return nothing
    o.truncLen > 0 && (last = o.truncLen)
    low = findfirst(q -> q <= o.truncQ, r.quality)
    low !== nothing && (last = min(last, low-1))
    # truncLen is measured before trimLeft, including after quality trimming.
    o.truncLen > 0 && last < o.truncLen && return nothing
    first = o.trimLeft + 1
    o.minLen <= last-first+1 <= o.maxLen || return nothing
    seq = r.sequence[first:last]
    count(==('N'), seq) <= o.maxN || return nothing
    q = @view r.quality[first:last]
    expectedErrors(q) <= o.maxEE || return nothing
    FastqRecord(r.id, seq, q)
end
function check_paths(inputs, outputs)
    ins = realpath.(inputs)
    outs = [joinpath(realpath(dirname(abspath(p))), basename(p)) for p in outputs]
    length(unique(outs)) == length(outs) || throw(ArgumentError("Outputs must differ"))
    any(p -> p in ins, outs) && throw(ArgumentError("Output would overwrite input"))
    any(isfile, outputs) && throw(ArgumentError("Output already exists"))
end
function filterAndTrim(input::AbstractString, output::AbstractString; options=FilterOptions())
    check_options(options)
    check_paths([input], [output])
    reader = FastqReader(input)
    n = kept = 0
    try
        out = open_sequence(output, "w")
        try
            for r in reader
                n += 1
                f = filtered(r, options)
                f === nothing && continue
                write_record(out, f)
                kept += 1
            end
        finally
            close(out)
        end
    catch
        isfile(output) && rm(output)
        rethrow()
    finally
        close(reader)
    end
    FilterStats(n, kept)
end
pair_id(id) = replace(first(split(id)), r"/[12]$" => "")
function filterAndTrim(forward::AbstractString, reverse::AbstractString, outF::AbstractString, outR::AbstractString;
                       optionsF=FilterOptions(), optionsR=FilterOptions())
    check_options(optionsF); check_options(optionsR)
    check_paths([forward, reverse], [outF, outR])
    resources = IO[]
    readers = FastqReader[]
    n = kept = 0
    created = String[]
    try
        push!(readers, FastqReader(forward)); push!(readers, FastqReader(reverse))
        push!(resources, open_sequence(outF, "w")); push!(created, outF)
        push!(resources, open_sequence(outR, "w")); push!(created, outR)
        while true
            f = iterate(readers[1]); r = iterate(readers[2])
            f === nothing && r === nothing && break
            (f === nothing || r === nothing) && throw(ArgumentError("Paired FASTQ lengths differ"))
            pair_id(f[1].id) == pair_id(r[1].id) || throw(ArgumentError("Paired read IDs differ"))
            n += 1
            ff = filtered(f[1], optionsF); rr = filtered(r[1], optionsR)
            (ff === nothing || rr === nothing) && continue
            write_record(resources[1], ff); write_record(resources[2], rr)
            kept += 1
        end
    catch
        foreach(close, resources); empty!(resources)
        foreach(p -> isfile(p) && rm(p), created)
        rethrow()
    finally
        foreach(close, resources); foreach(close, readers)
    end
    FilterStats(n, kept)
end
mutable struct DerepBuilder
    index::Dict{String,Int}
    sequences::Vector{String}
    counts::Vector{Int}
    quality_sums::Vector{Vector{Float64}}
    mapping::Vector{Int}
    keepmap::Bool
end
DerepBuilder(;keepmap=true) = DerepBuilder(Dict{String,Int}(), String[], Int[],
                                           Vector{Float64}[], Int[], keepmap)

function add_read!(builder::DerepBuilder, read::FastqRecord)
    seq = validated_dna(read.sequence)
    i = get(builder.index, seq, 0)
    if i == 0
        push!(builder.sequences, seq)
        push!(builder.counts, 0)
        push!(builder.quality_sums, zeros(length(read.quality)))
        i = length(builder.sequences)
        builder.index[seq] = i
    end
    builder.counts[i] += 1
    builder.quality_sums[i] .+= read.quality
    builder.keepmap && push!(builder.mapping, i)
    builder
end

function finish(builder::DerepBuilder)
    qualities = [s ./ builder.counts[i] for (i,s) in enumerate(builder.quality_sums)]
    Derep(builder.sequences, builder.counts, qualities, builder.mapping)
end

function dereplicate(records; keepmap=true)
    builder = DerepBuilder(;keepmap)
    foreach(read -> add_read!(builder,read), records)
    finish(builder)
end

"Dereplicate a whole-read prefix containing at least `target_bases` bases."
function dereplicate_prefix(path::AbstractString, target_bases::Integer)
    target_bases > 0 || throw(ArgumentError("target_bases must be positive"))
    builder = DerepBuilder(;keepmap=false)
    reads = bases = 0
    reader = FastqReader(path)
    try
        for read in reader
            add_read!(builder,read)
            reads += 1
            bases += length(read.quality)
            bases >= target_bases && break
        end
    finally
        close(reader)
    end
    finish(builder), reads, bases
end
function derepFastq(path::AbstractString; keepmap=true)
    reader = FastqReader(path)
    try
        dereplicate(reader; keepmap)
    finally
        close(reader)
    end
end
function reverseComplement(s::AbstractString)
    bytes = codeunits(validated_dna(s; iupac=true))
    source = codeunits("ACGTRYSWKMBDHVN"); target = codeunits("TGCAYRSWMKVHDBN")
    String([target[findfirst(==(b),source)] for b in reverse(bytes)])
end
reverseComplement(s::BioSequences.NucSeq)=BioSequences.reverse_complement(s)
