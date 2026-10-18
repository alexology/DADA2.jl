"Reusable sparse taxonomy index; create with `trainTaxonomy`."
struct TaxonomyModel
    k::Int
    taxonomy::Vector{Vector{String}}
    logtotals::Vector{Float64}
    postings::Dict{UInt32,Vector{Pair{Int,Float64}}}
end

"Train once and reuse the sparse classifier across samples or batches."
function trainTaxonomy(reference::AbstractString; k=8)
    1 <= k <= 10 || throw(ArgumentError("Invalid taxonomy k"))
    fasta = readFasta(reference)
    isempty(fasta.ids) && throw(ArgumentError("Empty taxonomy reference"))
    labels = [filter(!isempty,strip.(split(id,';'))) for id in fasta.ids]
    any(isempty,labels) && throw(ArgumentError("Empty taxonomy label"))
    # Group identical full taxonomic paths so scores do not depend on record order.
    paths = sort!(unique(join.(labels,";")))
    index = Dict(p => i for (i,p) in enumerate(paths))
    vocab = 1 << (2k)
    counts = [Dict{UInt32,Int}() for _ in paths]
    totals = fill(Float64(vocab),length(paths))
    for (s,label) in zip(fasta.sequences,labels)
        col = index[join(label,";")]
        for word in unique(kmers(s,k))
            counts[col][word] = get(counts[col],word,0)+1
            totals[col] += 1
        end
    end
    # Sparse inverted index: store only observed (word, taxon) combinations,
    # avoiding a 4^k × taxa dense training matrix.
    postings = Dict{UInt32,Vector{Pair{Int,Float64}}}()
    for col in eachindex(counts), (word,n) in counts[col]
        push!(get!(postings,word,Pair{Int,Float64}[]),col=>log1p(n))
    end
    logtotals = log.(totals)
    taxonomy = [String.(split(p,';')) for p in paths]
    TaxonomyModel(k,taxonomy,logtotals,postings)
end

"K-mer naive Bayes with bootstrap support and conservative common ancestry."
function assignTaxonomy(sequences::AbstractVector{<:DNAInput}, model::TaxonomyModel;
                        minBoot=50, nBoot=100, seed=1, tryRC=false)
    0 <= minBoot <= 100 && nBoot >= 1 || throw(ArgumentError("Invalid taxonomy options"))
    k = model.k; taxonomy = model.taxonomy
    logtotals = model.logtotals; postings = model.postings
    ranks = maximum(length,taxonomy)
    seq = [validated_dna(s; iupac=true) for s in sequences]
    taxa = Matrix{Union{Missing,String}}(missing,length(seq),ranks)
    boot = zeros(Int,length(seq),ranks)
    rng = MersenneTwister(seed)
    function classify(words)
        scores = -length(words) .* logtotals
        for w in words
            for (col,logcount) in get(postings,w,Pair{Int,Float64}[])
                scores[col] += logcount
            end
        end
        best = maximum(scores)
        winners = findall(x -> isapprox(x,best; atol=1e-10,rtol=0),scores)
        common = String[]
        for rank in 1:ranks
            all(c -> length(taxonomy[c]) >= rank,winners) || break
            values = [taxonomy[c][rank] for c in winners]
            length(unique(values)) == 1 || break
            push!(common,String(values[1]))
        end
        common, maximum(scores)
    end
    for (i,s) in enumerate(seq)
        words = unique(kmers(s,k)); isempty(words) && continue
        prediction, score = classify(words)
        if tryRC
            rcwords = unique(kmers(reverseComplement(s),k))
            rcpred, rcscore = classify(rcwords)
            if rcscore > score
                words = rcwords; prediction = rcpred
            end
        end
        supports = zeros(Int,length(prediction))
        for _ in 1:nBoot
            sampled = rand(rng,words,max(1,length(words)÷8))
            p, _ = classify(sampled)
            for rank in eachindex(prediction)
                rank <= length(p) && p[1:rank] == prediction[1:rank] && (supports[rank] += 1)
            end
        end
        for rank in eachindex(prediction)
            boot[i,rank] = round(Int,100supports[rank]/nBoot)
            if boot[i,rank] >= minBoot && (rank == 1 || !ismissing(taxa[i,rank-1]))
                taxa[i,rank] = prediction[rank]
            end
        end
    end
    TaxonomyResult(seq,taxa,boot)
end
assignTaxonomy(sequences::AbstractVector{<:DNAInput}, reference::AbstractString; k=8,kwargs...) =
    assignTaxonomy(sequences,trainTaxonomy(reference;k);kwargs...)

"Reusable seed index for exact species matches to full reference sequences."
struct SpeciesModel
    references::Vector{String}
    species::Vector{Set{String}}
    postings::Dict{UInt32,Vector{Int}}
    k::Int
end
function trainSpecies(reference::AbstractString; k=8)
    1 <= k <= 15 || throw(ArgumentError("Invalid species k"))
    fasta = readFasta(reference); lookup = Dict{String,Set{String}}()
    for (id,s) in zip(fasta.ids,fasta.sequences)
        fields = split(id)
        length(fields) >= 3 || throw(ArgumentError("Species headers must be: accession Genus species"))
        name = join(fields[2:3]," ")
        push!(get!(lookup,s,Set{String}()),name)
    end
    references = sort!(collect(keys(lookup)))
    postings = Dict{UInt32,Vector{Int}}()
    for (i,s) in enumerate(references), word in unique(kmers(s,k))
        push!(get!(postings,word,Int[]),i)
    end
    SpeciesModel(references,[lookup[s] for s in references],postings,k)
end

function species_matches(query,model::SpeciesModel)
    words = unique(kmers(query,model.k))
    candidates = eachindex(model.references)
    if !isempty(words)
        seed = words[argmin([length(get(model.postings,w,Int[])) for w in words])]
        candidates = get(model.postings,seed,Int[])
    end
    matches = Set{String}()
    for i in candidates
        occursin(query,model.references[i]) && union!(matches,model.species[i])
    end
    matches
end

"Exact full or fragment species matching; ambiguous species are missing by default."
function addSpecies(sequences::AbstractVector{<:DNAInput}, model::SpeciesModel; allowMultiple=false, tryRC=true)
    result = Vector{Union{Missing,String}}(missing,length(sequences))
    for (i,s) in enumerate(validated_dna.(sequences))
        isempty(s) && throw(ArgumentError("Empty species query"))
        matches = species_matches(s,model)
        tryRC && union!(matches,species_matches(reverseComplement(s),model))
        if length(matches) == 1 || (allowMultiple && !isempty(matches))
            result[i] = join(sort!(collect(matches)),"/")
        end
    end
    result
end
addSpecies(sequences::AbstractVector{<:DNAInput},reference::AbstractString;kwargs...) =
    addSpecies(sequences,trainSpecies(reference);kwargs...)
