# Outputs of the DADA2 1.36.0 C++ core, written by validation/core_reference.R.
const CORE = joinpath(@__DIR__, "fixtures", "core")
const CORE_OPTIONS = Dict("OMEGA_C"=>:omegaC, "GREEDY"=>:greedy, "GAPLESS"=>:gapless,
    "VECTORIZED_ALIGNMENT"=>:vectorizedAlignment, "HOMOPOLYMER_GAP_PENALTY"=>:homopolymerGap,
    "BAND_SIZE"=>:band, "USE_KMERS"=>:useKmers, "DETECT_SINGLETONS"=>:detectSingletons,
    "MIN_FOLD"=>:minFold, "MIN_HAMMING"=>:minHamming, "MIN_ABUNDANCE"=>:minAbundance,
    "OMEGA_A"=>:omegaA, "KDIST_CUTOFF"=>:kdistCutoff, "MAX_CLUST"=>:maxClusters)

core_value(s) = s == "TRUE" ? true : s == "FALSE" ? false :
    something(tryparse(Int, s), parse(Float64, s))
gapped(v) = String([b == 0 ? '-' : "ACGT"[b] for b in v])
read_tsv(path) = [split(line, '\t') for line in readlines(path)]

function core_case(dir)
    rows = read_tsv(joinpath(dir, "uniques.tsv"))
    derep = Derep(String.(getindex.(rows, 1)), parse.(Int, getindex.(rows, 2)),
                  [Float64.(codeunits(r[3])) .- 33 for r in rows])
    err = ErrorModel(permutedims(hcat([parse.(Float64, r) for r in read_tsv(joinpath(dir, "err.tsv"))]...)))
    kwargs = Dict{Symbol,Any}()
    for line in readlines(joinpath(dir, "options.txt"))
        key, value = split(line, '=')
        kwargs[CORE_OPTIONS[key]] = core_value(value)
    end
    priors = joinpath(dir, "priors.txt")
    isfile(priors) && (kwargs[:priors] = readlines(priors))
    denoised = read_tsv(joinpath(dir, "denoised.tsv"))[2:end]
    sequences = String.(getindex.(denoised, 1))
    expected = (asvs=Dict(zip(sequences, parse.(Int, getindex.(denoised, 2)))),
                map=[line == "NA" ? "" : sequences[parse(Int, line)] for line in readlines(joinpath(dir, "map.txt"))],
                trans=permutedims(hcat([parse.(Int, r) for r in read_tsv(joinpath(dir, "trans.tsv"))]...)))
    derep, err, kwargs, expected
end

@testset "Alignment identical to the upstream C++ aligners" begin
    rows = read_tsv(joinpath(CORE, "nwalign.tsv"))[2:end]
    @test length(rows) == 400
    for (s1, s2, band, endsfree, vec, homo, al0, al1) in rows
        a = DADA2.encode_dna(s1); b = DADA2.encode_dna(s2)
        band = parse(Int, band); free = endsfree == "TRUE"
        al = vec == "TRUE" ? DADA2.nwalign_vectorized(a, b, 5, -4, -8, free ? 0 : -8, band) :
            DADA2.nwalign_endsfree(a, b, 5, -4, -8, band; homo_gap=parse(Int, homo), endsfree=free)
        @test (gapped(al.al0), gapped(al.al1)) == (al0, al1)
    end
    gl = DADA2.nwalign_gapless(DADA2.encode_dna("ACGTA"), DADA2.encode_dna("ACG"))
    @test (gapped(gl.al0), gapped(gl.al1)) == ("ACGTA", "ACG--")
    sub = DADA2.al2subs(DADA2.PairAlignment(DADA2.encode_dna("ACGT"), UInt8[1, 0, 4, 4]))
    @test sub.map == [1, 0, 2, 3] && sub.pos == [3] && (sub.nt0, sub.nt1) == ([3], [4])
end

@testset "Partitions identical to the upstream C++ core" begin
    cases = filter(isdir, readdir(CORE; join=true))
    @test length(cases) == 14
    for dir in cases, threaded in (false, true)
        derep, err, kwargs, expected = core_case(dir)
        result = dada(derep, err; threaded, kwargs...)
        @test Dict(zip(result.sequences, result.counts)) == expected.asvs
        @test [i == 0 ? "" : result.sequences[i] for i in result.assignments] == expected.map
        @test result.transitions == expected.trans
        @test sum(result.counts) + result.rejected == sum(derep.counts)
    end
end

@testset "Core helpers" begin
    @test DADA2.calc_pA(1, 0.5, true) ≈ 1 - exp(-0.5)
    @test DADA2.calc_pA(2, 0.1, false) ≈ (1 - exp(-0.1)*1.1)/(1 - exp(-0.1))
    @test DADA2.calc_pA(2, 1e-9, false) ≈ 5e-10 rtol=1e-6  # Taylor normalisation
    # The ordered k-mer count equals the shared count only for gapless pairs.
    a = DADA2.kmer_words(DADA2.encode_dna("ACGTACGTTT")); b = DADA2.kmer_words(DADA2.encode_dna("ACGTACGATT"))
    @test DADA2.kord_dot(a, b) == DADA2.kmer_dot(sort(a), sort(b))
    c = DADA2.kmer_words(DADA2.encode_dna("CGTACGTTTA"))
    @test DADA2.kord_dot(a, c) < DADA2.kmer_dot(sort(a), sort(c))
    short = Derep(["ACGTA", "ACGTT"], [5, 1], [fill(30.0, 5), fill(30.0, 5)])
    @test_throws ArgumentError dada(short)
    # A cache built with other options is reset, not reused.
    d = Derep(["ACGTACGTACGT", "ACGTACGAACGT"], [50, 3], [fill(30.0, 12), fill(30.0, 12)])
    cache = DADA2.InferenceCache(d)
    first = dada(d; _cache=cache)
    @test dada(d; _cache=cache, band=-1).assignments == first.assignments
    @test cache.options.band == -1
end

@testset "Error learning identical to R learnErrors" begin
    root = joinpath(@__DIR__, "fixtures", "learn")
    files = [joinpath(root, "sample$i.fastq.gz") for i in 1:3]
    matrix(path) = permutedims(hcat([parse.(Float64, r) for r in read_tsv(path)]...))
    for name in ("all", "two")
        nbases, rounds = parse.(Int, readlines(joinpath(root, "$(name)_info.txt")))
        for threaded in (false, true)
            learned = learnErrors(files; nbases, threaded)
            @test learned.transitions == round.(Int, matrix(joinpath(root, "$(name)_trans.tsv")))
            @test learned.errors.rates ≈ matrix(joinpath(root, "$(name)_err.tsv")) rtol=1e-12
            @test learned.iterations == rounds && learned.converged
            @test learned.quality_mode == :mean
        end
    end
    derep = derepFastq.(files)
    @test learnErrors(derep).errors.rates == learnErrors(files).errors.rates
    learned = learnErrors(derep)
    @test [r.sequences for r in dada(derep, learned)] == [r.sequences for r in dada(derep, learned.errors)]
    @test_throws ArgumentError learnErrors(files; maxClusters=2)
end

@testset "Review regressions DD5-DD8" begin
    # DD5: Int16 scores overflow upstream from about 6,500 identical bases.
    rng = MersenneTwister(5)
    for len in (1000, 6600, 7000)
        s = String(rand(rng, UInt8.(['A','C','G','T']), len))
        d = Derep([s], [2], [fill(30.0, len)])
        @test sum(dada(d).transitions) == 2len
        @test dada(d).transitions == dada(d; vectorizedAlignment=false).transitions
        al = DADA2.nwalign_vectorized(DADA2.encode_dna(s), DADA2.encode_dna(s), 5, -4, -8, 0, 16)
        @test al.al0 == al.al1 == DADA2.encode_dna(s)
    end
    # DD6: every new comparison reaches the cache, also with several tasks.
    seqs = unique([String(rand(rng, UInt8.(['A','C','G','T']), 6)) for _ in 1:200])[1:64]
    d = Derep(seqs, fill(2, 64), [fill(30.0, 6) for _ in 1:64])
    complete = map(1:200) do _
        cache = DADA2.InferenceCache(d)
        dada(d; _cache=cache, threaded=true, maxClusters=1)
        length(cache.compare[cache.order[1]]) == 64
    end
    @test all(complete)
    # DD7: headers that readFasta would reject are refused before writing.
    mktempdir() do dir
        for id in ("", "   ")
            path = joinpath(dir, "empty.fasta")
            @test_throws ArgumentError writeFasta(path, ["ACGTACGT"]; ids=[id])
            @test !isfile(path)
        end
    end
    # DD8: the extended-precision sum ignores the caller's BigFloat precision.
    reference = DADA2.long_double_sum([0.1, 0.2, 0.3])
    @test setprecision(() -> DADA2.long_double_sum([0.1, 0.2, 0.3]), BigFloat, 32) == reference
    @test setprecision(() -> DADA2.long_double_sum([0.1, 0.2, 0.3]), BigFloat, 512) == reference
end

@testset "Fast paths equal the reference computations" begin
    rng = MersenneTwister(21)
    bases = UInt8.(['A','C','G','T'])
    for _ in 1:300
        # Long repeats push some k-mer counts above the 255 kept in the dense bytes.
        s = vcat(rand(rng, bases, rand(rng, 6:300)), repeat(rand(rng, bases, rand(rng, 1:3)), rand(rng, 0:400)))
        t = vcat(rand(rng, bases, rand(rng, 6:300)), repeat(rand(rng, bases, rand(rng, 1:3)), rand(rng, 0:400)))
        a = DADA2.RawSeq(DADA2.encode_dna(String(s)), UInt8[], true)
        b = DADA2.RawSeq(DADA2.encode_dna(String(t)), UInt8[], true)
        wa = DADA2.kmer_words(a.seq); wb = DADA2.kmer_words(b.seq)
        @test DADA2.kmer_dot(a, b) == DADA2.kmer_dot(sort(wa), sort(wb))
        @test DADA2.kord_dot(a.seq, b.seq) == DADA2.kord_dot(wa, wb)
        x = rand(rng, 0x01:0x04, rand(rng, 20:250)); y = copy(x)
        for _ in 1:rand(rng, 0:6)
            y[rand(rng, eachindex(y))] = rand(rng, 0x01:0x04)
        end
        rand(rng) < 0.4 && deleteat!(y, rand(rng, 2:length(y)-1))
        for (p, q) in ((x, y), (y, x)), band in (16, -1)
            full = DADA2.CompactSub(DADA2.nwalign_vectorized(p, q, 5, -4, -8, 0, band))
            fast = DADA2.nwalign_vectorized(p, q, 5, -4, -8, 0, band; compact=Val(true))
            @test (full.pos, full.trans) == (fast.pos, fast.trans)
        end
    end
end

@testset "Review regressions DD9-DD11" begin
    root = joinpath(@__DIR__, "fixtures", "learn")
    files = [joinpath(root, "sample$i.fastq.gz") for i in 1:3]
    dereps = derepFastq.(files)
    matrix(path) = permutedims(hcat([parse.(Float64, r) for r in read_tsv(path)]...))
    # DD10: dereplicated inputs use R's budget of whole samples.
    nbases = parse(Int, first(readlines(joinpath(root, "two_info.txt"))))
    budget = learnErrors(dereps; nbases)
    @test budget.transitions == round.(Int, matrix(joinpath(root, "two_trans.tsv")))
    @test budget.errors.rates ≈ matrix(joinpath(root, "two_err.tsv")) rtol=1e-12
    @test budget.transitions == learnErrors(files; nbases).transitions
    @test learnErrors(dereps[1]; nbases=1).transitions == learnErrors(dereps[1:1]).transitions
    @test_throws ArgumentError learnErrors(dereps; nbases=0)
    # DD11: quality sums are widened before they can overflow.
    builder = DADA2.DerepBuilder(; keepmap=false)
    DADA2.add_read!(builder, FastqRecord("a", "ACGTACGT", fill(40, 8)))
    DADA2.add_read!(builder, FastqRecord("b", "ACGTACGT", fill(40, 8)))
    limit = Int(typemax(UInt32) ÷ 93)
    builder.counts[1] = limit                       # state reached after `limit` reads
    builder.quality_sums[1] = fill(UInt32(40 * limit), 8)
    DADA2.add_read!(builder, FastqRecord("c", "ACGTACGT", fill(40, 8)))
    @test builder.quality_sums[1] isa Vector{UInt64}
    @test builder.quality_sums[1] == fill(UInt64(40) * (limit + 1), 8)
    @test DADA2.finish(builder).qualities[1] == fill(40.0, 8)
    # DD9: serial preparation gives the same uniques as the parallel one.
    d = dereps[1]
    serial = DADA2.InferenceCache(d; threaded=false); parallel = DADA2.InferenceCache(d; threaded=true)
    @test all(i -> (serial.raws[i].seq, serial.raws[i].qual, collect(serial.raws[i].kmer)) ==
                   (parallel.raws[i].seq, parallel.raws[i].qual, collect(parallel.raws[i].kmer)),
              eachindex(serial.raws))
    @test dada(d; threaded=false).assignments == dada(d; threaded=true).assignments
end
