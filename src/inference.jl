function phredErrors(; maxQ=93)
    0 <= maxQ <= 93 || throw(ArgumentError("maxQ must be in 0:93"))
    rates = zeros(16, maxQ+1)
    for q in 0:maxQ, from in 1:4, to in 1:4
        e = min(0.75, 10.0^(-q/10))
        rates[4(from-1)+to,q+1] = from == to ? 1-e : e/3
    end
    ErrorModel(rates)
end

"P(Poisson(expected) ≥ reads | X > 0), evaluated in log space."
function abundance_logp(reads::Int, expected::Float64; prior=false)
    reads >= 1 && isfinite(expected) && expected >= 0 || throw(ArgumentError("Invalid Poisson arguments"))
    expected == 0 && return (reads == 1 && !prior ? 0.0 : -Inf)
    lp = logccdf(Poisson(expected), reads-1)
    prior || (lp -= log(-expm1(-expected)))
    min(0.0, lp)
end
abundancePvalue(reads::Integer, expected::Real; prior=false) = exp(abundance_logp(Int(reads), Float64(expected); prior))

"""
    dada(derep, err; omegaA=1e-40, omegaP=1e-4, omegaC=1e-40, priors=String[], ...)

Divisive partitioning of the unique sequences, ported from the DADA2 1.36.0
C++ core with the defaults of `setDadaOpt`: ends-free banded alignment
(`band=16`, `match=5`, `mismatch=-4`, `gap=-8`), k-mer screen
(`kdistCutoff=0.42`), gapless shortcut, greedy locking, abundance p-values
with Bonferroni correction and the final `omegaC` correction. Reads whose
uniques fail the final test are left unassigned (`assignments` zero).
`homopolymerGap` enables homopolymer gap penalties (non-vectorized aligner);
`minFold`, `minHamming` and `minAbundance` restrict which sequences can start
a new partition. `converged` is false when the maximum number of shuffles or
`maxClusters` stopped the algorithm early.
"""
function dada(d::Derep, err::ErrorModel=phredErrors(); omegaA=1e-40, omegaP=1e-4, omegaC=1e-40,
              priors=String[], detectSingletons=false, band=16, useKmers=true,
              kdistCutoff=0.42, maxClusters=0, match=5, mismatch=-4, gap=-8,
              homopolymerGap=nothing, vectorizedAlignment=true, gapless=true, greedy=true,
              minFold=1, minHamming=1, minAbundance=1, threaded=true,
              cacheAlignments=true, _cache=nothing)
    0 < omegaA <= 1 && 0 < omegaP <= 1 && 0 <= omegaC <= 1 && 0 <= kdistCutoff <= 1 &&
        maxClusters >= 0 && match > 0 && mismatch <= 0 && minFold >= 0 &&
        all(x -> abs(x) <= 1000, (match, mismatch, gap)) ||
        throw(ArgumentError("Invalid inference options"))
    # Positive penalties are negated, as in R.
    gap = -abs(Int(gap))
    homo = homopolymerGap === nothing ? gap : -abs(Int(homopolymerGap))
    abs(homo) <= 1000 || throw(ArgumentError("Invalid inference options"))
    vectorized = vectorizedAlignment && homo == gap && band != 0
    maxq = size(err.rates,2)-1
    all(q -> all(x -> round(Int, x, RoundNearestTiesAway) <= maxq, q), d.qualities) ||
        throw(ArgumentError("Error model does not cover read qualities"))
    n = length(d.sequences)
    n == 0 && return DadaResult(String[], Int[], Int[], d, zeros(Int,size(err.rates)), err, true, 0, 0,
                                empty_clustering())
    all(s -> KMER_SIZE < ncodeunits(s) < 9999, d.sequences) ||
        throw(ArgumentError("Sequences must be longer than $KMER_SIZE and shorter than 9999 bases"))
    cache = _cache === nothing ? InferenceCache(d; useKmers, cacheAlignments) : _cache
    options = AlignOptions(match, mismatch, gap, homo, Int(band), useKmers, kdistCutoff,
                           vectorized, gapless)
    check_cache!(cache, d, options, cacheAlignments)
    prior_set = Set(validated_dna.(priors))
    order = cache.order
    b = Partition(cache.raws, d.counts[order], [d.sequences[i] in prior_set for i in order],
                  Float64(omegaA), Float64(omegaP))
    o = DadaOptions(options, err.rates, greedy, detectSingletons, Float64(minFold),
                    Int(minHamming), Int(minAbundance), threaded)
    converged = run_dada!(b, cache, o, Int(maxClusters))

    # Final alignments with each center, then the unconditional correction test.
    members = [(c.center, r) for c in b.clusters for r in c.raws]
    subs = Vector{Sub}(undef, n)
    foreach_scratch(length(members), threaded) do scratch, k
        center, r = members[k]
        subs[r] = final_sub(b, cache, o, center, r, scratch)
    end
    for c in b.clusters, r in c.raws
        if r == c.center
            b.p[r] = 1.0
        else
            b.p[r] = calc_pA(b.reads[r], b.comp[r].lambda * c.reads, true)
            b.p[r] < omegaC && (b.correct[r] = false)
        end
    end

    nclust = length(b.clusters)
    representative = zeros(Int, nclust); totals = zeros(Int, nclust)
    assignments = zeros(Int, n); rejected = 0
    trans = zeros(Int, size(err.rates))
    for (i, c) in enumerate(b.clusters)
        max_reads = 0
        center = b.raws[c.center]
        for r in c.raws
            if b.reads[r] > max_reads
                representative[i] = r; max_reads = b.reads[r]
            end
            if !b.correct[r]
                rejected += b.reads[r]
                continue
            end
            totals[i] += b.reads[r]; assignments[order[r]] = i
            raw = b.raws[r]; map = subs[r].map
            for pos0 in eachindex(center.seq)
                pos1 = map[pos0]
                pos1 == 0 && continue
                t = 4*(center.seq[pos0]-1) + raw.seq[pos1]
                trans[t, raw.qual[pos1]+1] += b.reads[r]
            end
        end
    end
    sequences = [d.sequences[order[r]] for r in representative]
    table = clustering_table(b, cache, o, subs, sequences, totals)
    ranked = sortperm(1:nclust; by=i -> (-totals[i], sequences[i]))
    inverse = invperm(ranked)
    remapped = [i == 0 ? 0 : inverse[i] for i in assignments]
    DadaResult(sequences[ranked], totals[ranked], remapped, d, trans, err, converged, b.nalign, rejected,
               select_clustering(table, ranked, inverse))
end

"Upstream `b_make_clustering_df`, partitions in order of birth."
function clustering_table(b::Partition, cache::InferenceCache, o::DadaOptions, subs, sequences, totals)
    nclust = length(b.clusters)
    n0 = zeros(Int, nclust); n1 = zeros(Int, nclust); nunq = zeros(Int, nclust)
    for (i, c) in enumerate(b.clusters), r in c.raws
        b.correct[r] || continue
        nunq[i] += 1
        k = nsubs(subs[r])
        k == 0 && (n0[i] += b.reads[r])
        k == 1 && (n1[i] += b.reads[r])
    end
    # Post-hoc p-value: reads expected at each centre from the other partitions.
    center_of = Dict(c.center => i for (i, c) in enumerate(b.clusters))
    expected = zeros(nclust)
    for (i, c) in enumerate(b.clusters), comp in c.comps
        j = get(center_of, Int(comp.index), 0)
        j != 0 && j != i && (expected[j] += comp.lambda * c.reads)
    end
    pval = [calc_pA(b.reads[c.center], expected[i], true) for (i, c) in enumerate(b.clusters)]
    birth_from = Union{Missing,Int}[missing]; birth_pval = Union{Missing,Float64}[missing]
    birth_fold = Union{Missing,Float64}[missing]; birth_ham = Union{Missing,Int}[missing]
    birth_qave = Union{Missing,Float64}[missing]
    for c in b.clusters[2:end]
        parent = b.clusters[c.birth_comp.i].center
        sub = sub_new(b.raws[parent], b.raws[c.center], o.align, o.align.use_kmers, 1.0)
        qave = 0.0
        if sub !== nothing
            q = b.raws[c.center].qual
            qave = sum(Float64(q[sub.map[p]]) for p in sub.pos; init=0.0) / nsubs(sub)
        end
        push!(birth_from, c.birth_from); push!(birth_pval, c.birth_pval)
        push!(birth_fold, c.birth_fold); push!(birth_ham, c.birth_comp.hamming)
        push!(birth_qave, qave)
    end
    ClusteringTable((sequences, totals, n0, n1, nunq, pval, birth_from, birth_pval, birth_fold,
                     birth_ham, birth_qave))
end

"Query position → center position (0 for an insertion), from processing indices."
function center_mapping(cache::InferenceCache, c::Int, r::Int)
    cache.options === nothing && throw(ArgumentError("Inference cache was not used"))
    entry = cache.final[r]
    sub = entry !== nothing && entry[1] == c ? entry[2] :
          sub_new(cache.raws[c], cache.raws[r], cache.options, false, 1.0)
    mapping = zeros(Int, length(cache.raws[r].seq))
    for (pos0, pos1) in enumerate(sub.map)
        pos1 == 0 || (mapping[pos1] = pos0)
    end
    mapping
end

dada(d::Derep, err::AbstractMatrix; kwargs...) = dada(d, ErrorModel(err); kwargs...)
dada(path::AbstractString, err=phredErrors(); kwargs...) = dada(derepFastq(path), err; kwargs...)
# Always finish sibling tasks before returning or propagating their first error.
function wait_all(tasks)
    first_error = nothing
    for task in tasks
        try
            wait(task)
        catch e
            first_error === nothing && (first_error = e)
        end
    end
    first_error === nothing || throw(first_error)
    nothing
end

"""
    dada(samples, err; pool=false, priors=String[], pseudoPrevalence=2, pseudoAbundance=Inf, ...)

Infer each sample. `pool=true` denoises the pooled uniques once (R's
`combineDereps2`) and returns each sample's share of the pooled partitions;
`pool=:pseudo` first infers each sample, then infers again with the ASVs found
in at least `pseudoPrevalence` samples or `pseudoAbundance` reads as priors.
As in R, that second pass uses the error model re-estimated by `loessErrfun`
from the transitions of the first pass, not `err`.
"""
function dada(samples::AbstractVector{Derep}, err::ErrorModel=phredErrors(); pool=false,
              priors=String[], pseudoPrevalence=2, pseudoAbundance=Inf, kwargs...)
    pool in (false, true, :pseudo) || throw(ArgumentError("pool must be false, true or :pseudo"))
    length(samples) <= 1 && (pool = false)
    if pool === true
        haskey(kwargs, :_caches) && throw(ArgumentError("Pooled inference does not use sample caches"))
        pooled, index = pool_dereps(samples)
        return split_pooled(dada(pooled, err; priors, kwargs...), samples, index)
    end
    results = dada_each(samples, fill(err, length(samples)); priors, kwargs...)
    pool === :pseudo || return results
    pseudoPrevalence >= 0 && pseudoAbundance >= 0 || throw(ArgumentError("Invalid pseudo-pooling options"))
    prevalence = Dict{String,Int}(); abundance = Dict{String,Int}()
    for r in results, (sequence, count) in zip(r.sequences, r.counts)
        count > 0 || continue
        prevalence[sequence] = get(prevalence, sequence, 0) + 1
        abundance[sequence] = get(abundance, sequence, 0) + count
    end
    pseudo = [s for s in keys(prevalence) if prevalence[s] >= pseudoPrevalence || abundance[s] >= pseudoAbundance]
    # As upstream, the second pass uses the error model re-estimated with
    # loessErrfun from the transitions of the first pass.
    err2 = r_loess_errors(accumulate_transitions(results), nothing)
    errs = [widen_errors(err2, ceil(Int, maximum(maximum, d.qualities; init=0.0))) for d in samples]
    dada_each(samples, errs; priors=vcat(collect(priors), sort!(pseudo)), kwargs...)
end

"R `accumulateTrans`: transition counts summed, columns padded to the widest."
function accumulate_transitions(results)
    trans = zeros(Int, 16, maximum(r -> size(r.transitions, 2), results))
    for r in results
        trans[:, 1:size(r.transitions, 2)] .+= r.transitions
    end
    trans
end

"R's extension of an error model to quality `q`: the last column is repeated."
widen_errors(err::ErrorModel, q::Int) = size(err.rates, 2) >= q + 1 ? err :
    ErrorModel(hcat(err.rates, repeat(err.rates[:, end:end], 1, q + 1 - size(err.rates, 2))))

"Upstream `combineDereps2`: summed counts, weighted mean qualities, by decreasing abundance."
function pool_dereps(samples::AbstractVector{Derep})
    # Each sample in its processing order (decreasing abundance, as R derep objects).
    sequences = String[]; position = Dict{String,Int}()
    for d in samples, i in sortperm(1:length(d.sequences); by=i -> (-d.counts[i], i))
        haskey(position, d.sequences[i]) && continue
        push!(sequences, d.sequences[i]); position[d.sequences[i]] = length(sequences)
    end
    counts = zeros(Int, length(sequences))
    sums = [zeros(ncodeunits(s)) for s in sequences]
    for d in samples, i in sortperm(1:length(d.sequences); by=i -> (-d.counts[i], i))
        j = position[d.sequences[i]]
        counts[j] += d.counts[i]
        sums[j] .+= d.qualities[i] .* d.counts[i]
    end
    qualities = [sums[j] ./ counts[j] for j in eachindex(sums)]
    order = sortperm(1:length(sequences); by=j -> (-counts[j], j))
    pooled = Derep(sequences[order], counts[order], qualities[order])
    pooled, Dict(s => k for (k, s) in enumerate(pooled.sequences))
end

"R's expansion of a pooled result: the partitions holding each sample's reads."
function split_pooled(pooled::DadaResult, samples, index)
    map(samples) do d
        counts = zeros(Int, length(pooled.sequences)); rejected = 0
        for (i, s) in enumerate(d.sequences)
            c = pooled.assignments[index[s]]
            c == 0 ? (rejected += d.counts[i]) : (counts[c] += d.counts[i])
        end
        present = unique(c for c in (pooled.assignments[index[s]] for s in d.sequences) if c != 0)
        kept = sort(present; by=c -> (-counts[c], pooled.sequences[c]))
        renumber = zeros(Int, length(counts)); renumber[kept] = eachindex(kept)
        assignments = [(c = pooled.assignments[index[s]]; c == 0 ? 0 : renumber[c]) for s in d.sequences]
        table = select_clustering(pooled.clustering, kept, renumber)
        table.abundance .= counts[kept]
        DadaResult(pooled.sequences[kept], counts[kept], assignments, d, pooled.transitions,
                   pooled.errors, pooled.converged, pooled.comparisons, rejected, table)
    end
end

"Independent inference of each sample with its own error model."
function dada_each(samples::AbstractVector{Derep}, errs::AbstractVector{ErrorModel}; threaded=true,
                   _caches=nothing, kwargs...)
    length(errs) == length(samples) || throw(ArgumentError("One error model per sample is required"))
    results = Vector{DadaResult}(undef,length(samples))
    caches = _caches === nothing ? fill(nothing,length(samples)) : _caches
    length(caches) == length(samples) || throw(ArgumentError("Invalid inference caches"))
    if threaded && Threads.nthreads() > 1 && length(samples) > 1
        # Independent tasks let the scheduler redistribute uneven samples.
        tasks = map(eachindex(samples)) do i
            Threads.@spawn dada(samples[i],errs[i];threaded=false,_cache=caches[i],kwargs...)
        end
        wait_all(tasks)
        for i in eachindex(tasks)
            results[i] = fetch(tasks[i])
        end
    else
        # With one sample there is no outer work to distribute. Preserve
        # threading for the comparisons inside `dada` instead of disabling it.
        inner_threaded = threaded && length(samples) == 1
        for i in eachindex(samples)
            results[i] = dada(samples[i],errs[i]; threaded=inner_threaded,_cache=caches[i],kwargs...)
        end
    end
    results
end

const SMOOTHERS = (:loess, :direct, :linear)

"Previous local linear estimator, retained for explicit scientific comparisons."
function linear_fitted_errors(trans::Matrix{Int}, fallback::ErrorModel; span=0.75)
    rates = copy(fallback.rates); nq = size(trans,2)
    for from in 1:4
        rows = 4(from-1)+1:4from
        totals = vec(sum(@view(trans[rows,:]); dims=1))
        observed = findall(>(0),totals)
        isempty(observed) && continue
        neighbors = max(2,ceil(Int,span*length(observed)))
        for to in 1:4
            from == to && continue
            row = 4(from-1)+to
            y = [log10((trans[row,q]+1)/totals[q]) for q in observed]
            for q in 1:nq
                distances = abs.(observed .- q)
                radius = sort(distances)[min(neighbors,length(distances))] + 1
                sw = sx = sy = sxx = sxy = 0.0
                for j in eachindex(observed)
                    x = observed[j]-q
                    w = totals[observed[j]] * max(0.0,1-(abs(x)/radius)^3)^3
                    sw += w; sx += w*x; sy += w*y[j]; sxx += w*x*x; sxy += w*x*y[j]
                end
                denom = sw*sxx-sx*sx
                pred = denom > eps()*sw*max(sxx,1) ? (sy*sxx-sx*sxy)/denom : sy/sw
                rates[row,q] = clamp(10.0^pred,1e-7,0.25)
            end
        end
        for q in 1:nq
            rates[4(from-1)+from,q] = 1-sum(rates[4(from-1)+to,q] for to in 1:4 if to != from)
        end
    end
    ErrorModel(rates)
end

"""
Smoothed error rates from transition counts. `:loess` is DADA2's `loessErrfun`
(R's interpolated LOESS surface); `:direct` fits the local quadratic at every
quality (R's `surface="direct"`); `:linear` is the earlier local linear fit.
"""
function fitted_errors(trans::Matrix{Int}, fallback::ErrorModel; span=0.75, smoother=:loess)
    size(trans) == size(fallback.rates) && all(>=(0), trans) ||
        throw(ArgumentError("Invalid transition counts"))
    isfinite(span) && 0 < span <= 1 || throw(ArgumentError("Invalid smoothing span"))
    smoother in SMOOTHERS || throw(ArgumentError("smoother must be :loess, :direct or :linear"))
    smoother == :linear && return linear_fitted_errors(trans, fallback; span)
    smoother == :loess &&
        return r_loess_errors(trans, fallback; span, sparse=direct_fitted_errors(trans, fallback; span))
    direct_fitted_errors(trans, fallback; span)
end

# R `loess(..., surface="direct")`; supports of at most three qualities are
# interpolated linearly, where a local quadratic is not identifiable.
function direct_fitted_errors(trans::Matrix{Int}, fallback::ErrorModel; span=0.75)
    rates = copy(fallback.rates)
    nq = size(trans, 2)
    for from in 1:4
        totals = vec(sum(@view(trans[4(from-1)+1:4from, :]); dims=1))
        observed = findall(>(0), totals)
        isempty(observed) && continue # No evidence for this source: preserve initial rates.
        alternatives = [4(from-1)+to for to in 1:4 if to != from]
        responses = [log10((trans[row, q]+1)/totals[q]) for q in observed, row in alternatives]
        # R's direct LOESS uses floor(span*n) nearest points, tricube weights,
        # weighted least squares, and degree two. Small supports need lower degree.
        neighbors = clamp(floor(Int, span*length(observed)), 1, length(observed))
        for q in first(observed):last(observed)
            if length(observed) <= 3
                # Sparse/binned support cannot sustain a local quadratic fit.
                left = searchsortedlast(observed,q)
                right = min(left+1,length(observed))
                fraction = right == left ? 0.0 : (q-observed[left])/(observed[right]-observed[left])
                for (j,row) in enumerate(alternatives)
                    prediction = (1-fraction)*responses[left,j]+fraction*responses[right,j]
                    rates[row,q] = clamp(10.0^prediction,1e-7,0.25)
                end
                continue
            end
            distances = abs.(observed .- q)
            radius = Float64(sort(distances)[neighbors])
            # At very small spans the nearest point(s) provide the only evidence.
            weights = radius == 0 ? Float64.(distances .== 0) :
                [max(0.0, 1-(d/radius)^3)^3 for d in distances]
            if !any(>(0), weights)
                weights .= distances .== minimum(distances)
            end
            weights .*= totals[observed]
            weights ./= maximum(weights)
            active = findall(>(0), weights)
            degree = min(2, length(active)-1)
            scale = max(radius, 1.0)
            x = (observed[active] .- q) ./ scale
            w = sqrt.(weights[active])
            design = [w[j]*x[j]^power for j in eachindex(active), power in 0:degree]
            # Reuse the same small factorization for all three destination bases.
            coefficients = pinv(design) * (responses[active, :] .* w)
            for (j, row) in enumerate(alternatives)
                rates[row, q] = clamp(10.0^coefficients[1,j], 1e-7, 0.25)
            end
        end
        for row in alternatives
            rates[row, 1:first(observed)-1] .= rates[row, first(observed)]
            rates[row, last(observed)+1:nq] .= rates[row, last(observed)]
        end
        for q in 1:nq
            rates[4(from-1)+from, q] = 1-sum(rates[row,q] for row in alternatives)
        end
    end
    ErrorModel(rates)
end

struct ErrorLearningResult
    errors::ErrorModel
    history::Vector{Float64}
    converged::Bool
    iterations::Int
    transitions::Matrix{Int}
    inference_converged::Bool
    quality_mode::Symbol
    training_reads::Int
    training_bases::Int
end

struct ErrorTrainingSource
    path::String
    reads::Int
    bases::Int
end

function raw_transitions_sample!(trans, d, result, source::ErrorTrainingSource, cache::InferenceCache;
                                 band=16)
    source.reads == 0 && return trans
    index = Dict(s=>i for (i,s) in enumerate(d.sequences))
    reference_indices = zeros(Int,length(d.sequences))
    for i in eachindex(d.sequences)
        c = result.assignments[i]
        c == 0 && continue
        reference_indices[i] = index[result.sequences[c]]
    end
    seen = zeros(Int,length(d.sequences))
    position = invperm(cache.order)
    mappings = Vector{Union{Nothing,Vector{Int}}}(nothing,length(d.sequences))
    reader = FastqReader(source.path)
    nread = 0
    try
        for read in reader
            nread += 1
            i = get(index,read.sequence,0)
            i > 0 || throw(ArgumentError("FASTQ changed during error learning"))
            seen[i] += 1
            c = result.assignments[i]
            if c != 0
                ref = result.sequences[c]
                center = reference_indices[i]
                if center != i && mappings[i] === nothing
                    mappings[i] = center_mapping(cache,position[center],position[i])
                end
                mapping = center == i ? nothing : mappings[i]
                for pos in eachindex(read.quality)
                    mapped = mapping === nothing ? pos : mapping[pos]
                    mapped == 0 && continue
                    from = baseindex(codeunit(ref,mapped)); to = baseindex(codeunit(read.sequence,pos))
                    q = Int(read.quality[pos])+1
                    q <= size(trans,2) || throw(ArgumentError("Error model does not cover raw read qualities"))
                    trans[4(from-1)+to,q] += 1
                end
            end
            nread == source.reads && break
        end
    finally
        close(reader)
    end
    nread == source.reads && seen == d.counts ||
        throw(ArgumentError("FASTQ changed during error learning"))
    trans
end

function raw_transitions!(trans, samples, results, sources, caches; band=16, threaded=true)
    length(samples) == length(results) == length(sources) == length(caches) ||
        throw(ArgumentError("Training sources and samples differ"))
    partial = [zeros(Int,size(trans)) for _ in eachindex(samples)]
    if threaded && Threads.nthreads() > 1 && length(samples) > 1
        tasks = map(eachindex(samples)) do i
            Threads.@spawn raw_transitions_sample!(partial[i],samples[i],results[i],sources[i],caches[i];band)
        end
        wait_all(tasks)
    else
        for i in eachindex(samples)
            raw_transitions_sample!(partial[i],samples[i],results[i],sources[i],caches[i];band)
        end
    end
    fill!(trans,0)
    for counts in partial
        trans .+= counts
    end
    trans
end

# R `dada(..., selfConsist=TRUE)` as called by `learnErrors`. With
# `rawsources`, transitions are recounted at the original read qualities.
function learn_errors_upstream(samples::AbstractVector{Derep}; maxIter=10, span=0.75, omegaC=0.0,
                               threaded=true, cacheAlignments=true, useKmers=true,
                               rawsources=nothing, rawqmax=nothing, kwargs...)
    !isempty(samples) && all(s -> !isempty(s.counts), samples) ||
        throw(ArgumentError("Cannot learn from empty samples"))
    maxIter >= 1 && isfinite(span) && 0 < span <= 1 || throw(ArgumentError("Invalid learning options"))
    haskey(kwargs, :maxClusters) && throw(ArgumentError("maxClusters is set by error learning"))
    caches = [InferenceCache(sample; useKmers, cacheAlignments) for sample in samples]
    # R: qmax = ceiling(max(quals)); models are widened to cover it. Raw
    # counting needs the highest original quality as well.
    qmax = [ceil(Int, maximum(maximum, s.qualities)) for s in samples]
    rawqmax === nothing || (qmax = max.(qmax, rawqmax))
    widen(err, q) = size(err.rates, 2) >= q + 1 ? err :
        ErrorModel(hcat(err.rates, repeat(err.rates[:, end:end], 1, q + 1 - size(err.rates, 2))),
                   Val(:unchecked))  # the first model does not sum to one
    function round_trans(errs, maxClusters)
        results = dada_each(samples, errs; omegaC, threaded, _caches=caches, maxClusters, useKmers,
                            cacheAlignments, kwargs...)
        converged = all(r -> r.converged, results)
        rawsources === nothing && return accumulate_transitions(results), converged
        trans = zeros(Int, 16, maximum(e -> size(e.rates, 2), errs))
        raw_transitions!(trans, samples, results, rawsources, caches; threaded), converged
    end
    # Initialisation: maximal error rates and a single partition per sample.
    initial = [ErrorModel(ones(16, max(41, q + 1)), Val(:unchecked)) for q in qmax]
    trans, _ = round_trans(initial, 1)
    err = r_loess_errors(trans, nothing; span)
    rates = copy(err.rates)
    rates[[1, 6, 11, 16], :] .= 1.0
    err = ErrorModel(rates, Val(:unchecked))
    previous = ErrorModel[]; history = Float64[]
    converged = false; inference_converged = false
    for _ in 1:maxIter
        push!(previous, err)
        trans, inference_converged = round_trans([widen(err, q) for q in qmax], 0)
        updated = r_loess_errors(trans, nothing; span)
        push!(history, size(updated.rates) == size(err.rates) ?
            maximum(abs(log(max(updated.rates[i],1e-300)) - log(max(err.rates[i],1e-300)))
                    for i in eachindex(err.rates)) : Inf)
        err = updated
        # R stops when a model repeats exactly, or after MAX_CONSIST rounds.
        if any(e -> e.rates == err.rates, previous)
            converged = true; break
        end
    end
    training_reads = sum(sum(s.counts) for s in samples)
    training_bases = sum(sum(length(s.sequences[i])*s.counts[i] for i in eachindex(s.sequences))
                         for s in samples)
    ErrorLearningResult(err, history, converged, length(history), trans, inference_converged,
                        rawsources === nothing ? :mean : :raw, training_reads, training_bases)
end

"""
R `learnErrors` sampling: whole files in order until the bases exceed `nbases`.
Returns the dereplicated files, their read sources and highest original quality.
"""
function upstream_training_data(paths::AbstractVector{<:AbstractString}, nbases::Integer)
    !isempty(paths) || throw(ArgumentError("At least one FASTQ path is required"))
    nbases > 0 || throw(ArgumentError("nbases must be positive"))
    samples = Derep[]; sources = ErrorTrainingSource[]; rawqmax = Int[]; total = 0
    for path in paths
        builder = DerepBuilder(; keepmap=false)
        reads = bases = qmax = 0
        reader = FastqReader(path)
        try
            for read in reader
                add_read!(builder, read)
                reads += 1; bases += length(read.quality)
                qmax = max(qmax, Int(maximum(read.quality; init=0x00)))
            end
        finally
            close(reader)
        end
        push!(samples, finish(builder)); push!(sources, ErrorTrainingSource(String(path), reads, bases))
        push!(rawqmax, qmax)
        total += bases
        total > nbases && break
    end
    samples, sources, rawqmax
end

"""
    learnErrors(samples_or_paths; nbases=100_000_000, maxIter=10, qualityMode=:mean, ...)

Learn one error model shared by the samples, as R `learnErrors`: whole FASTQ
files in order until more than `nbases` bases, a first round with maximal
error rates and one partition per sample, self-transitions set to one in that
first model, then up to `maxIter=10` rounds of inference (`omegaC=0`) and
DADA2's `loessErrfun` that stop when a model repeats exactly. With the same
input it gives R's error model, to rounding. Other keywords go to `dada`.

Transitions are counted at the mean qualities of the dereplicated reads, as
R. For FASTQ inputs, `qualityMode=:raw` counts them at the original read
qualities instead (Julia extension; same procedure otherwise, files reread
at each round). Use `errorDiagnostics(result)` to inspect training support;
convergence does not demonstrate external calibration or biological accuracy.
`cacheAlignments=false` reduces persistent memory use for unusually large
numbers of unique sequences, at the cost of repeating alignments each round.
"""
function learnErrors end

# Options of the learner removed in favour of R's procedure.
const REMOVED_LEARNER_OPTIONS = (:initial, :tolerance, :smoother, :upstream, :rawsources, :rawqmax)

function check_learner_options(kwargs)
    for key in REMOVED_LEARNER_OPTIONS
        haskey(kwargs, key) && throw(ArgumentError("learnErrors follows R learnErrors; `$key` is not supported"))
    end
end

function learnErrors(samples::AbstractVector{Derep}; qualityMode=:mean, kwargs...)
    qualityMode == :mean || throw(ArgumentError("Raw quality learning requires FASTQ paths"))
    check_learner_options(kwargs)
    learn_errors_upstream(samples; kwargs...)
end
learnErrors(d::Derep; kwargs...) = learnErrors([d]; kwargs...)

function learnErrors(paths::AbstractVector{<:AbstractString}; qualityMode=:mean,
                     nbases::Integer=100_000_000, kwargs...)
    qualityMode in (:mean, :raw) || throw(ArgumentError("qualityMode must be :mean or :raw"))
    check_learner_options(kwargs)
    samples, sources, rawqmax = upstream_training_data(paths, nbases)
    qualityMode == :mean && return learn_errors_upstream(samples; kwargs...)
    learn_errors_upstream(samples; rawsources=sources, rawqmax, kwargs...)
end
learnErrors(path::AbstractString; kwargs...) = learnErrors([path]; kwargs...)

"Training support by source base and quality; convergence is not external calibration."
function errorDiagnostics(learned::ErrorLearningResult)
    nq = size(learned.transitions,2)
    exposures = zeros(Int,4,nq); observed = zeros(Int,4,nq); expected = zeros(4,nq)
    for from in 1:4, q in 1:nq
        rows = 4(from-1)+1:4from
        exposures[from,q] = sum(@view learned.transitions[rows,q])
        observed[from,q] = exposures[from,q]-learned.transitions[4(from-1)+from,q]
        expected[from,q] = exposures[from,q] *
            sum(learned.errors.rates[4(from-1)+to,q] for to in 1:4 if to != from)
    end
    (qualities=collect(0:nq-1), exposures=exposures, observed_substitutions=observed,
     expected_substitutions=expected, supported=exposures .> 0,
     unobserved_sources=["ACGT"[b] for b in 1:4 if all(iszero,@view exposures[b,:])],
     converged=learned.converged, inference_converged=learned.inference_converged,
     quality_mode=learned.quality_mode, training_reads=learned.training_reads,
     training_bases=learned.training_bases)
end
dada(d::Derep, learned::ErrorLearningResult; kwargs...) = dada(d,learned.errors; kwargs...)
dada(samples::AbstractVector{Derep}, learned::ErrorLearningResult; kwargs...) =
    dada(samples, learned.errors; kwargs...)
