"Streaming four-line FASTQ reader (FASTX.jl); caller owns and closes the reader."
mutable struct FastqReader
    reader::FASTX.FASTQ.Reader
    record::Int
end
function open_sequence(path::AbstractString, mode::AbstractString; gzip=endswith(lowercase(path), ".gz"))
    raw = open(path, mode)
    if gzip
        try
            return mode == "r" ? GzipDecompressorStream(raw) : GzipCompressorStream(raw)
        catch
            close(raw)
            rethrow()
        end
    end
    raw
end
FastqReader(path::AbstractString) = FastqReader(FASTX.FASTQReader(open_sequence(path, "r"); copy=false), 0)
Base.close(r::FastqReader) = close(r.reader)
Base.IteratorSize(::Type{FastqReader}) = Base.SizeUnknown()
Base.eltype(::Type{FastqReader}) = FastqRecord
function Base.iterate(r::FastqReader, state=nothing)
    next = try
        iterate(r.reader)
    catch e
        e isa ErrorException || e isa EOFError || rethrow()
        throw(ArgumentError("Invalid FASTQ after record $(r.record): $(sprint(showerror, e))"))
    end
    next === nothing && return nothing
    rec = first(next)
    r.record += 1
    qual = [Int(c)-33 for c in codeunits(FASTX.quality(rec))]
    (FastqRecord(FASTX.description(rec), FASTX.sequence(String, rec), qual), nothing)
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
# Same order as R `fastqFilter`: maxLen on the raw read, trimLeft, truncQ on
# the bases kept after trimLeft, then truncLen (measured on the raw read; R
# ignores it when it does not exceed trimLeft) and minLen.
function filtered(r::FastqRecord, o::FilterOptions)
    last = length(r.quality)
    last <= o.maxLen || return nothing
    first = o.trimLeft + 1
    last >= first || return nothing
    low = findnext(q -> q <= o.truncQ, r.quality, first)
    low !== nothing && (last = low-1)
    if o.truncLen >= first
        last < o.truncLen && return nothing
        last = o.truncLen
    end
    last-first+1 >= o.minLen || return nothing
    seq = r.sequence[first:last]
    count(==('N'), seq) <= o.maxN || return nothing
    q = @view r.quality[first:last]
    expectedErrors(q) <= o.maxEE || return nothing
    FastqRecord(r.id, seq, q)
end
# On a case-insensitive file system (Windows, macOS by default) names that
# differ only in letter case are the same file; probed on the directory itself.
function case_insensitive(dir)
    up, low = uppercase(dir), lowercase(dir)
    up != low && ispath(up) && ispath(low) && samefile(up, low)
end
function path_key(p)
    dir = realpath(dirname(abspath(p)))
    key = joinpath(dir, basename(p))
    case_insensitive(dir) ? lowercase(key) : key
end
function check_paths(inputs, outputs)
    for p in inputs
        isfile(p) || throw(ArgumentError("Input file not found: $p"))
    end
    for p in outputs
        dir = dirname(abspath(p))
        isdir(dir) || throw(ArgumentError("Output directory does not exist: $dir"))
    end
    ins = path_key.(realpath.(inputs))
    outs = path_key.(outputs)
    length(unique(outs)) == length(outs) || throw(ArgumentError("Outputs must differ"))
    any(p -> p in ins, outs) && throw(ArgumentError("Output would overwrite input"))
    any(isfile, outputs) && throw(ArgumentError("Output already exists"))
end
function filterAndTrim(input::AbstractString, output::AbstractString; options=FilterOptions())
    check_options(options)
    check_paths([input], [output])
    reader = FastqReader(input)
    n = kept = 0
    tmp = nothing
    try
        tmp = temp_output(output)
        out = open_sequence(tmp, "w"; gzip=endswith(lowercase(output), ".gz"))
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
        publish_outputs([tmp], [output])
    catch
        tmp === nothing || discard(tmp)
        rethrow()
    finally
        close(reader)
    end
    FilterStats(n, kept)
end
# Outputs are written to temporary files in the destination directory and
# renamed into place only after the whole run succeeded, so a failure never
# leaves a partial file under the requested name.
function temp_output(path)
    tmp, io = mktemp(dirname(abspath(path)); cleanup=false)
    close(io)
    tmp
end
# Renames are not a transaction: if a later one fails, the outputs already
# published by this call are removed again. A single rename either publishes
# the whole file or leaves the destination untouched; `mv` would fall back to
# copying and then fail, leaving the output in place despite the error.
function publish_outputs(tmps, outputs)
    done = String[]
    try
        for (t, o) in zip(tmps, outputs)
            # Never overwrite a file that appeared after `check_paths`.
            ispath(o) && throw(ArgumentError("Output already exists: $o"))
            # Julia 1.10's Filesystem.rename falls back to copying on failure.
            err = ccall(:jl_fs_rename, Int32, (Cstring, Cstring), t, o)
            err < 0 && Base.uv_error("rename($(repr(t)), $(repr(o)))", err)
            push!(done, o)
        end
    catch
        foreach(discard, done)
        rethrow()
    end
end
# Best-effort removal during cleanup: a failure here must not hide the
# original error.
function discard(path)
    try
        rm(path; force=true)
    catch e
        @warn "Could not remove $path" exception=e
    end
end
pair_id(id) = replace(first(split(id)), r"/[12]$" => "")
function filterAndTrim(forward::AbstractString, reverse::AbstractString, outF::AbstractString, outR::AbstractString;
                       optionsF=FilterOptions(), optionsR=FilterOptions())
    check_options(optionsF); check_options(optionsR)
    check_paths([forward, reverse], [outF, outR])
    resources = IO[]
    readers = FastqReader[]
    n = kept = 0
    tmps = String[]
    try
        push!(readers, FastqReader(forward)); push!(readers, FastqReader(reverse))
        for o in (outF, outR)
            push!(tmps, temp_output(o))
            push!(resources, open_sequence(tmps[end], "w"; gzip=endswith(lowercase(o), ".gz")))
        end
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
        foreach(close, resources); empty!(resources)
        publish_outputs(tmps, [outF, outR])
    catch
        foreach(close, resources); empty!(resources)
        foreach(discard, tmps)
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
