"Merge denoised paired reads using ungapped suffix/prefix overlaps."
function mergePairs(forward::DadaResult, reverse::DadaResult; minOverlap=12, maxMismatch=0)
    minOverlap >= 1 && maxMismatch >= 0 || throw(ArgumentError("Invalid merge options"))
    fm = forward.derep.map; rm = reverse.derep.map
    isempty(fm) && throw(ArgumentError("mergePairs requires keepmap=true"))
    length(fm) == length(rm) || throw(ArgumentError("Paired read counts differ"))
    pairs = Dict{Tuple{Int,Int},Int}()
    for i in eachindex(fm)
        key = (forward.assignments[fm[i]], reverse.assignments[rm[i]])
        pairs[key] = get(pairs,key,0)+1
    end
    sequences = Dict{String,Int}(); accepted = rejected = 0
    for ((fi,ri),count) in sort!(collect(pairs); by=first)
        if fi == 0 || ri == 0
            rejected += count
            continue
        end
        f = forward.sequences[fi]; r = reverseComplement(reverse.sequences[ri])
        best = nothing; bestmismatch = typemax(Int)
        for overlap in min(ncodeunits(f),ncodeunits(r)):-1:minOverlap
            mismatches = 0
            for j in 1:overlap
                mismatches += codeunit(f,ncodeunits(f)-overlap+j) != codeunit(r,j)
                mismatches > maxMismatch && break
            end
            if mismatches <= maxMismatch && mismatches < bestmismatch
                best = f * r[overlap+1:end]; bestmismatch = mismatches
                mismatches == 0 && break
            end
        end
        if best === nothing
            rejected += count
        else
            accepted += count; sequences[best] = get(sequences,best,0)+count
        end
    end
    order = sort!(collect(keys(sequences)); by=s -> (-sequences[s],s))
    MergeResult(order, [sequences[s] for s in order], accepted, rejected)
end

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

"Exact two-parent chimera test; parents must be strictly more abundant."
function isBimeraDenovo(sequence::DNAInput, parents::AbstractVector{<:DNAInput}; minSegment=8)
    minSegment >= 1 || throw(ArgumentError("minSegment must be positive"))
    query = validated_dna(sequence); n = ncodeunits(query)
    prefix = Int[]; suffix = Int[]
    for parent in parents
        p = validated_dna(parent)
        p == query && continue
        ncodeunits(p) == n || continue
        left = right = 0
        while left < n && codeunit(query,left+1) == codeunit(p,left+1)
            left += 1
        end
        while right < n && codeunit(query,n-right) == codeunit(p,n-right)
            right += 1
        end
        push!(prefix,left); push!(suffix,right)
    end
    # Only the two longest suffixes matter: one may share the prefix's parent.
    best = second = 0
    for j in eachindex(suffix)
        if best == 0 || suffix[j] > suffix[best]
            second = best; best = j
        elseif second == 0 || suffix[j] > suffix[second]
            second = j
        end
    end
    for i in eachindex(prefix)
        j = best == i ? second : best
        j == 0 && continue
        max(minSegment,n-suffix[j]) <= min(n-minSegment,prefix[i]) && return true
    end
    false
end
function removeBimeraDenovo(table::SequenceTable; minFoldParentOverAbundance=1.0, minSegment=8,
                          method=:consensus, minSampleFraction=0.9)
    minFoldParentOverAbundance >= 1 && 0 <= minSampleFraction <= 1 && minSegment >= 1 || throw(ArgumentError("Invalid chimera options"))
    method in (:consensus,:pooled) || throw(ArgumentError("method must be :consensus or :pooled"))
    keep = trues(length(table.sequences))
    totals = vec(sum(table.counts; dims=1))
    for j in eachindex(keep)
        if method == :pooled
            parents = [table.sequences[k] for k in eachindex(keep) if k != j && totals[k] > minFoldParentOverAbundance*totals[j]]
            keep[j] = !isBimeraDenovo(table.sequences[j],parents; minSegment)
        else
            tested = votes = 0
            for sample in axes(table.counts,1)
                abundance = table.counts[sample,j]
                abundance == 0 && continue
                tested += 1
                parents = [table.sequences[k] for k in eachindex(keep) if k != j && table.counts[sample,k] > minFoldParentOverAbundance*abundance]
                votes += isBimeraDenovo(table.sequences[j],parents; minSegment)
            end
            keep[j] = !(votes > 0 && votes/tested >= minSampleFraction)
        end
    end
    SequenceTable(table.samples, table.sequences[keep], table.counts[:,keep])
end

function readFasta(path::AbstractString; bio=false)
    io = open_sequence(path,"r"); ids = String[]; sequences = String[]; buffer = IOBuffer()
    try
        for line in eachline(io)
            isempty(strip(line)) && continue
            if startswith(line,">")
                if !isempty(ids)
                    s = validated_dna(String(take!(buffer)); iupac=true)
                    isempty(s) && throw(ArgumentError("Empty FASTA sequence"))
                    push!(sequences,s)
                end
                id = strip(line[2:end])
                isempty(id) && throw(ArgumentError("Empty FASTA header"))
                push!(ids,id)
            else
                isempty(ids) && throw(ArgumentError("FASTA sequence before header"))
                write(buffer,strip(line))
            end
        end
        if !isempty(ids)
            s = validated_dna(String(take!(buffer)); iupac=true)
            isempty(s) && throw(ArgumentError("Empty FASTA sequence"))
            push!(sequences,s)
        end
    finally
        close(io)
    end
    (ids=ids,sequences=bio ? BioSequences.LongDNA{4}.(sequences) : sequences)
end
function writeFasta(path::AbstractString, sequences::AbstractVector{<:DNAInput}; ids=["ASV$i" for i in eachindex(sequences)])
    length(ids) == length(sequences) || throw(ArgumentError("Header count differs"))
    isfile(path) && throw(ArgumentError("Output already exists"))
    seq = [validated_dna(s; iupac=true) for s in sequences]
    any(isempty,seq) && throw(ArgumentError("Empty FASTA sequence"))
    all(id -> !occursin('\n',id) && !occursin('\r',id),ids) || throw(ArgumentError("Invalid FASTA header"))
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
