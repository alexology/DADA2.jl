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

function transition_logprob(ref, query, qualities, err, mapping)
    a = codeunits(ref); b = codeunits(query)
    lp = 0.0
    for pos in eachindex(b)
        from = mapping[pos] == 0 ? baseindex(b[pos]) : baseindex(a[mapping[pos]])
        to = baseindex(b[pos]); q = floor(Int, qualities[pos]+0.5) + 1
        lp += err.logrates[4(from-1)+to,q]
    end
    lp
end

const CachedMapping = Union{Nothing,Vector{Int}}
const CenterMappings = Vector{CachedMapping}

"Reusable, training-local sequence features and pairwise alignment mappings."
mutable struct InferenceCache
    spectra::Vector{Vector{UInt32}}
    mappings::Vector{Union{Nothing,CenterMappings}}
    band::Int
    cache_alignments::Bool
end

function InferenceCache(d::Derep; useKmers=true, band=16, cacheAlignments=true)
    spectra = useKmers ? kmers.(d.sequences) : Vector{UInt32}[]
    mappings = Vector{Union{Nothing,CenterMappings}}(nothing,length(d.sequences))
    InferenceCache(spectra,mappings,band,cacheAlignments)
end

function cached_mapping!(cache::InferenceCache, d::Derep, center::Int, query::Int)
    center == query && return nothing
    cache.cache_alignments ||
        return alignSequences(d.sequences[center],d.sequences[query];band=cache.band).mapping
    mappings = cache.mappings[center]
    if mappings === nothing
        mappings = CenterMappings(nothing,length(d.sequences))
        cache.mappings[center] = mappings
    end
    mapping = mappings[query]
    if mapping === nothing
        mapping = alignSequences(d.sequences[center],d.sequences[query];band=cache.band).mapping
        mappings[query] = mapping
    end
    mapping
end

"Infer variants by repeated abundance tests and likelihood-based repartitioning."
function dada(d::Derep, err::ErrorModel=phredErrors(); omegaA=1e-40, omegaP=1e-4, omegaC=1e-40,
              priors=String[], detectSingletons=false, band=16, useKmers=true,
              kdistCutoff=0.42, maxClusters=0, maxConsist=20, threaded=true,
              cacheAlignments=true, _cache=nothing)
    0 < omegaA <= 1 && 0 < omegaP <= 1 && 0 <= omegaC <= 1 && 0 <= kdistCutoff <= 1 &&
        band >= 0 && maxClusters >= 0 && maxConsist >= 1 || throw(ArgumentError("Invalid inference options"))
    maxq = size(err.rates,2)-1
    all(q -> all(x -> floor(Int,x+0.5) <= maxq, q), d.qualities) || throw(ArgumentError("Error model does not cover read qualities"))
    n = length(d.sequences)
    n == 0 && return DadaResult(String[], Int[], Int[], d, zeros(Int,size(err.rates)), err, true, 0, 0)
    prior_set = Set(validated_dna.(priors))
    centers = [sortperm(1:n; by=i -> (-d.counts[i], d.sequences[i]))[1]]
    iscenter = falses(n); iscenter[centers[1]] = true
    assignments = ones(Int,n); totals = [sum(d.counts)]
    cache = _cache === nothing ? InferenceCache(d;useKmers,band,cacheAlignments) : _cache
    cache.band == band && length(cache.mappings) == n &&
        cache.cache_alignments == cacheAlignments &&
        (!useKmers || length(cache.spectra) == n) || throw(ArgumentError("Invalid inference cache"))
    spectra = cache.spectra
    # One vector per discovered center; never allocate an n × n distance matrix.
    probabilities = Vector{Float64}[]
    comparison_count = 0
    function compare_center(center)
        out = fill(-Inf, n)
        # Allocate the center row before entering the threaded loop. Each
        # worker then owns a distinct query slot in that row.
        cache.cache_alignments && cache.mappings[center] === nothing &&
            (cache.mappings[center] = CenterMappings(nothing,n))
        function compare(i)
            if center == i
                out[i] = transition_logprob(d.sequences[i],d.sequences[i],d.qualities[i],err,eachindex(d.qualities[i]))
                return
            end
            if useKmers && center != i && kmer_distance(spectra[center], spectra[i]) > kdistCutoff
                return
            end
            mapping = cached_mapping!(cache,d,center,i)
            out[i] = transition_logprob(d.sequences[center], d.sequences[i], d.qualities[i], err, mapping)
        end
        if threaded && n >= 64 && Threads.nthreads() > 1
            Threads.@threads for i in 1:n
                compare(i)
            end
        else
            for i in 1:n
                compare(i)
            end
        end
        out
    end
    push!(probabilities, compare_center(centers[1])); comparison_count += n
    converged = true
    while true
        candidate = 0; bestp = Inf
        for i in 1:n
            iscenter[i] && continue
            prior = d.sequences[i] in prior_set
            d.counts[i] == 1 && !prior && !detectSingletons && continue
            c = assignments[i]
            lp = probabilities[c][i]
            expected = exp(lp + log(totals[c]))
            p = abundance_logp(d.counts[i], expected; prior=prior || detectSingletons)
            threshold = prior ? log(omegaP) : log(omegaA)-log(n)
            if p < threshold && (p < bestp || (p == bestp && (candidate == 0 || d.counts[i] > d.counts[candidate])))
                candidate = i; bestp = p
            end
        end
        candidate == 0 && break
        if maxClusters > 0 && length(centers) >= maxClusters
            converged = false; break
        end
        push!(centers, candidate); push!(probabilities, compare_center(candidate)); comparison_count += n
        iscenter[candidate] = true
        assignments[candidate] = length(centers)
        totals = zeros(Int,length(centers))
        for i in 1:n
            totals[assignments[i]] += d.counts[i]
        end
        stable = false
        for _ in 1:maxConsist
            next = copy(assignments)
            logtotals = log.(totals)
            for i in 1:n
                iscenter[i] && continue
                best = -Inf
                for c in eachindex(centers)
                    e = probabilities[c][i] + logtotals[c]
                    if e > best
                        best = e; next[i] = c
                    end
                end
            end
            if next == assignments
                stable = true; break
            end
            assignments = next; fill!(totals,0)
            for i in 1:n
                totals[assignments[i]] += d.counts[i]
            end
        end
        converged &= stable
        !stable && break
    end
    # Final unconditional abundance test: unsupported residual reads remain
    # unassigned (zero), as in upstream OMEGA_C error correction.
    rejected = 0
    output_totals = zeros(Int,length(centers))
    for i in 1:n
        c = assignments[i]
        expected = exp(probabilities[c][i] + log(totals[c]))
        if !iscenter[i] && omegaC > 0 && abundance_logp(d.counts[i],expected;prior=true) < log(omegaC)
            rejected += d.counts[i]; assignments[i] = 0
        else
            output_totals[c] += d.counts[i]
        end
    end
    totals = output_totals
    trans = zeros(Int, size(err.rates))
    for i in 1:n
        assignments[i] == 0 && continue
        ref = d.sequences[centers[assignments[i]]]
        center = centers[assignments[i]]
        mapping = iscenter[i] ? eachindex(d.qualities[i]) : cached_mapping!(cache,d,center,i)
        for pos in eachindex(mapping)
            mapping[pos] == 0 && continue
            from = baseindex(codeunit(ref, mapping[pos])); to = baseindex(codeunit(d.sequences[i],pos))
            q = floor(Int,d.qualities[i][pos]+0.5)+1
            trans[4(from-1)+to,q] += d.counts[i]
        end
    end
    order = sortperm(eachindex(centers); by=c -> (-totals[c], d.sequences[centers[c]]))
    inverse = invperm(order)
    remapped = [c == 0 ? 0 : inverse[c] for c in assignments]
    DadaResult(d.sequences[centers[order]], totals[order], remapped, d, trans, err, converged, comparison_count, rejected)
end
dada(d::Derep, err::AbstractMatrix; kwargs...) = dada(d, ErrorModel(err); kwargs...)
dada(path::AbstractString, err=phredErrors(); kwargs...) = dada(derepFastq(path), err; kwargs...)
function dada(samples::AbstractVector{Derep}, err::ErrorModel=phredErrors(); threaded=true,
              _caches=nothing, kwargs...)
    results = Vector{DadaResult}(undef,length(samples))
    caches = _caches === nothing ? fill(nothing,length(samples)) : _caches
    length(caches) == length(samples) || throw(ArgumentError("Invalid inference caches"))
    if threaded && Threads.nthreads() > 1 && length(samples) > 1
        # Independent tasks let the scheduler redistribute uneven samples.
        tasks = map(eachindex(samples)) do i
            Threads.@spawn dada(samples[i],err;threaded=false,_cache=caches[i],kwargs...)
        end
        for i in eachindex(tasks)
            results[i] = fetch(tasks[i])
        end
    else
        # With one sample there is no outer work to distribute. Preserve
        # threading for the comparisons inside `dada` instead of disabling it.
        inner_threaded = threaded && length(samples) == 1
        for i in eachindex(samples)
            results[i] = dada(samples[i],err; threaded=inner_threaded,_cache=caches[i],kwargs...)
        end
    end
    results
end

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

"Native direct quadratic LOESS; R's default uses an interpolated LOESS surface."
function fitted_errors(trans::Matrix{Int}, fallback::ErrorModel; span=0.75, smoother=:loess)
    size(trans) == size(fallback.rates) && all(>=(0), trans) ||
        throw(ArgumentError("Invalid transition counts"))
    isfinite(span) && 0 < span <= 1 || throw(ArgumentError("Invalid smoothing span"))
    smoother in (:loess, :linear) || throw(ArgumentError("smoother must be :loess or :linear"))
    smoother == :linear && return linear_fitted_errors(trans, fallback; span)
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
                mapping = center == i ? nothing : cached_mapping!(cache,d,center,i)
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
        for task in tasks
            fetch(task)
        end
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

function learn_errors(samples::AbstractVector{Derep}; initial=phredErrors(), maxIter=10, tolerance=1e-3,
                      span=0.75, smoother=:loess, rawsources=nothing, omegaC=0.0, kwargs...)
    !isempty(samples) && any(s -> !isempty(s.counts),samples) || throw(ArgumentError("Cannot learn from empty samples"))
    maxIter >= 1 && isfinite(tolerance) && tolerance > 0 && isfinite(span) && 0 < span <= 1 &&
        smoother in (:loess, :linear) || throw(ArgumentError("Invalid learning options"))
    err = initial; history = Float64[]; converged = false
    trans = zeros(Int,size(err.rates)); inference_converged = false
    use_kmers = get(kwargs,:useKmers,true)
    band = get(kwargs,:band,16)
    cache_alignments = get(kwargs,:cacheAlignments,true)
    caches = [InferenceCache(sample;useKmers=use_kmers,band,
        cacheAlignments=cache_alignments) for sample in samples]
    for _ in 1:maxIter
        results = dada(samples,err; omegaC,_caches=caches,kwargs...)
        fill!(trans,0)
        if rawsources === nothing
            for result in results
                trans .+= result.transitions
            end
        else
            raw_transitions!(trans,samples,results,rawsources,caches;band,
                             threaded=get(kwargs,:threaded,true))
        end
        inference_converged = all(r -> r.converged,results)
        updated = fitted_errors(trans,err; span, smoother)
        delta = maximum(abs(log(max(updated.rates[i],1e-300))-log(max(err.rates[i],1e-300))) for i in eachindex(err.rates))
        push!(history,delta); err = updated
        if delta < tolerance && inference_converged
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
    learnErrors(samples_or_paths; smoother=:loess, maxIter=10, ...)

Learn one shared substitution model while inferring each sample independently.
Native quadratic LOESS uses exposure weights, a +1 error pseudocount, and
constant rates outside observed qualities. `smoother=:linear` selects the
previous estimator for comparisons. Error learning uses `omegaC=0` by default
so every assigned read contributes to the transition counts.

FASTQ inputs default to `qualityMode=:raw`: use up to `nbases=100_000_000`
training bases and reread the selected prefixes each iteration to count
transitions at original qualities. `qualityMode=:mean` reproduces
dereplicated mean-quality counting. Derep inputs only support mean qualities.
Use `errorDiagnostics(result)` to inspect training support. Numerical
convergence does not demonstrate external calibration or biological accuracy.
`cacheAlignments=false` reduces persistent memory use for unusually large
numbers of unique sequences, at the cost of repeating alignments each iteration.
"""
function learnErrors(samples::AbstractVector{Derep}; qualityMode=:mean, kwargs...)
    qualityMode == :mean || throw(ArgumentError("Raw quality learning requires FASTQ paths"))
    learn_errors(samples;kwargs...)
end

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
learnErrors(d::Derep; kwargs...) = learnErrors([d]; kwargs...)
function error_training_data(paths::AbstractVector{<:AbstractString}, nbases::Integer)
    !isempty(paths) || throw(ArgumentError("At least one FASTQ path is required"))
    nbases > 0 || throw(ArgumentError("nbases must be positive"))
    samples = Derep[]
    sources = ErrorTrainingSource[]
    remaining = Int(nbases)
    for (i,path) in enumerate(paths)
        remaining == 0 && break
        target = cld(remaining,length(paths)-i+1)
        sample, reads, bases = dereplicate_prefix(path,target)
        reads == 0 && continue
        push!(samples,sample)
        push!(sources,ErrorTrainingSource(String(path),reads,bases))
        remaining = max(0,remaining-bases)
    end
    !isempty(samples) || throw(ArgumentError("Cannot learn from empty FASTQ files"))
    samples, sources
end

function learnErrors(paths::AbstractVector{<:AbstractString}; qualityMode=:raw,
                     nbases::Integer=100_000_000, kwargs...)
    qualityMode in (:raw,:mean) || throw(ArgumentError("qualityMode must be :raw or :mean"))
    samples, sources = error_training_data(paths,nbases)
    learn_errors(samples;rawsources=qualityMode == :raw ? sources : nothing,kwargs...)
end
learnErrors(path::AbstractString; kwargs...) = learnErrors([path]; kwargs...)
dada(d::Derep, learned::ErrorLearningResult; kwargs...) = dada(d,learned.errors; kwargs...)
