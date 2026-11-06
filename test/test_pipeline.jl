# R dada2 1.36.0 outputs written by validation/core_reference.R and
# validation/pipeline_reference.R: $clustering, pooling, mergePairs, chimeras.
const PIPELINE = joinpath(@__DIR__, "fixtures", "pipeline")

function uniques_derep(path; readmap=nothing)
    rows = read_tsv(path)
    map = readmap === nothing ? Int[] : parse.(Int, readlines(readmap))
    Derep(String.(getindex.(rows, 1)), parse.(Int, getindex.(rows, 2)),
          [Float64.(codeunits(r[3])) .- 33 for r in rows], map)
end
r_number(s) = s == "NA" ? missing : parse(Float64, s)
same_number(a, b) = (a === missing && b === missing) ||
    (a !== missing && b !== missing && (isnan(a) && isnan(b) || isapprox(a, b; rtol=1e-9)))

"Compare a clustering table with R's, rows matched by sequence."
function check_clustering(result::DadaResult, path; birth=true)
    rows = read_tsv(path)
    header = rows[1]; col(name) = findfirst(==(name), header)
    rseq = String.(getindex.(rows[2:end], col("sequence")))
    t = result.clustering
    @test sort(rseq) == sort(t.sequence)
    for row in rows[2:end]
        j = findfirst(==(row[col("sequence")]), t.sequence)
        j === nothing && continue
        @test t.abundance[j] == parse(Int, row[col("abundance")])
        @test (t.n0[j], t.n1[j], t.nunq[j]) ==
              Tuple(parse(Int, row[col(c)]) for c in ("n0", "n1", "nunq"))
        @test same_number(t.pval[j], r_number(row[col("pval")]))
        for c in ("birth_pval", "birth_fold", "birth_qave")
            @test same_number(getproperty(t, Symbol(c))[j], r_number(row[col(c)]))
        end
        ham = row[col("birth_ham")]
        @test t.birth_ham[j] === (ham == "NA" ? missing : parse(Int, ham))
        # Upstream leaves birth_from unset for births from priors (birth_pval >= omegaA).
        prior_birth = row[col("birth_pval")] != "NA" && parse(Float64, row[col("birth_pval")]) >= 1e-40
        if birth && !prior_birth
            from = row[col("birth_from")]
            parent = t.birth_from[j]
            @test (parent === missing ? "NA" : t.sequence[parent]) ==
                  (from == "NA" ? "NA" : rseq[parse(Int, from)])
        end
    end
end

function check_denoised(result::DadaResult, prefix)
    denoised = read_tsv(prefix * "_denoised.tsv")[2:end]
    sequences = String.(getindex.(denoised, 1))
    @test Dict(zip(result.sequences, result.counts)) == Dict(zip(sequences, parse.(Int, getindex.(denoised, 2))))
    expected = [l == "NA" ? "" : sequences[parse(Int, l)] for l in readlines(prefix * "_map.txt")]
    @test [i == 0 ? "" : result.sequences[i] for i in result.assignments] == expected
end

@testset "Clustering diagnostics identical to R" begin
    for dir in filter(isdir, readdir(CORE; join=true))
        derep, err, kwargs, _ = core_case(dir)
        check_clustering(dada(derep, err; kwargs...), joinpath(dir, "clustering.tsv"))
    end
end

@testset "Pooled and pseudo-pooled inference identical to R" begin
    err = ErrorModel(permutedims(hcat([parse.(Float64, r) for r in read_tsv(joinpath(PIPELINE, "err.tsv"))]...)))
    samples = [uniques_derep(joinpath(PIPELINE, "pool$(s)_uniques.tsv")) for s in 1:4]
    for (mode, pool) in (("pooled", true), ("pseudo", :pseudo)), threaded in (false, true)
        results = dada(samples, err; pool, threaded)
        for s in 1:4
            prefix = joinpath(PIPELINE, "$mode$s")
            check_denoised(results[s], prefix)
            # Pooled birth_from refers to the pooled partitions in R.
            check_clustering(results[s], prefix * "_clustering.tsv"; birth=mode == "pseudo")
            @test sum(results[s].counts) + results[s].rejected == sum(samples[s].counts)
        end
    end
    @test_throws ArgumentError dada(samples, err; pool=:unknown)
    @test dada(samples[1:1], err; pool=true)[1].sequences == dada(samples[1], err).sequences
end

@testset "mergePairs identical to R" begin
    err = ErrorModel(permutedims(hcat([parse.(Float64, r) for r in read_tsv(joinpath(PIPELINE, "err.tsv"))]...)))
    results = map(("F", "R")) do side
        d = uniques_derep(joinpath(PIPELINE, "merge$(side)_uniques.tsv");
                          readmap=joinpath(PIPELINE, "merge$(side)_readmap.txt"))
        r = dada(d, err)
        check_denoised(r, joinpath(PIPELINE, "merge$side"))
        check_clustering(r, joinpath(PIPELINE, "merge$(side)_clustering.tsv"))
        r
    end
    settings = Dict("default" => (;), "mismatch" => (maxMismatch=2, minOverlap=20),
                    "trim" => (trimOverhang=true,), "concat" => (justConcatenate=true,))
    for (name, options) in settings
        merged = mergePairs(results...; options...)
        p = merged.pairs
        julia_rows = [(p.sequence[k], p.abundance[k], results[1].sequences[p.forward[k]],
                       results[2].sequences[p.reverse[k]], p.nmatch[k], p.nmismatch[k], p.nindel[k],
                       p.prefer[k] === missing ? "NA" : string(p.prefer[k]), p.accept[k])
                      for k in eachindex(p.sequence)]
        r_rows = [(r[1], parse(Int, r[2]), r[3], r[4], parse(Int, r[5]), parse(Int, r[6]), parse(Int, r[7]),
                   r[8], r[9] == "TRUE") for r in read_tsv(joinpath(PIPELINE, "merge_$name.tsv"))[2:end]]
        @test julia_rows == r_rows
        @test merged.accepted == sum(p.abundance[p.accept]; init=0)
        @test merged.accepted + merged.rejected == length(results[1].derep.map)
    end
end

@testset "Chimera detection identical to R" begin
    rows = read_tsv(joinpath(PIPELINE, "chimera_table.tsv"))
    samples = String.(rows[1][2:end])
    sequences = String.(getindex.(rows[2:end], 1))
    counts = hcat([parse.(Int, r[2:end]) for r in rows[2:end]]...)
    table = SequenceTable(samples, sequences, counts)
    flags = read_tsv(joinpath(PIPELINE, "chimera_flags.tsv"))
    expected(name) = [r[findfirst(==(name), flags[1])] == "TRUE" for r in flags[2:end]]
    for threaded in (false, true)
        @test isBimeraDenovoTable(table; threaded) == expected("consensus")
        @test isBimeraDenovoTable(table; allowOneOff=true, threaded) == expected("consensus_oneoff")
        @test isBimeraDenovoTable(table; minSampleFraction=0.5, ignoreNNegatives=2,
            minFoldParentOverAbundance=1, minParentAbundance=4, maxShift=4, threaded) == expected("consensus_loose")
        @test isBimeraDenovo(table; threaded) == expected("pooled")
        @test isBimeraDenovo(table; allowOneOff=true, minOneOffParentDistance=2, threaded) == expected("pooled_oneoff")
    end
    @test removeBimeraDenovo(table).sequences == sequences[.!expected("consensus")]
    persample = read_tsv(joinpath(PIPELINE, "chimera_persample.tsv"))
    clean = removeBimeraDenovo(table; method=:per_sample)
    @test Dict(zip(clean.sequences, eachcol(clean.counts))) ==
          Dict(String(r[1]) => parse.(Int, r[2:end]) for r in persample[2:end])
end

@testset "PhiX screening identical to R" begin
    rows = read_tsv(joinpath(@__DIR__, "fixtures", "phix.tsv"))
    header = rows[1]; rows = rows[2:end]
    sequences = String.(getindex.(rows, 1))
    expected(name) = [r[findfirst(==(name), header)] == "TRUE" for r in rows]
    @test isPhiX(sequences) == expected("default")
    @test isPhiX(sequences; wordSize=12, minMatches=3) == expected("w12_m3")
    @test isPhiX(sequences; nonOverlapping=false, minMatches=3) == expected("overlapping")
    @test isPhiX(sequences; wordSize=20, minMatches=1) == expected("w20_m1")
    # filterAndTrim discards phiX reads unless rmPhix=false.
    mktempdir() do dir
        phix = DADA2.PHIX_GENOME[2001:2150]
        input = joinpath(dir, "in.fastq")
        write(input, "@phix\n$phix\n+\n$(repeat('I', 150))\n@other\n$(repeat("ACGTTGCA", 18))AC\n+\n$(repeat('I', 146))\n")
        @test filterAndTrim(input, joinpath(dir, "a.fastq")).reads_out == 1
        @test filterAndTrim(input, joinpath(dir, "b.fastq"); options=FilterOptions(rmPhix=false)).reads_out == 2
    end
end

@testset "Dereplication order as R derepFastq" begin
    reads = [FastqRecord("r$i", s, fill(30, 8)) for (i, s) in
             enumerate(["TTTTACGT", "AAAAACGT", "CCCCACGT", "AAAAACGT", "GGGGACGT", "CCCCACGT"])]
    # Decreasing abundance; ties alphabetical within the chunk where they first appear.
    d = dereplicate(reads)
    @test d.sequences == ["AAAAACGT", "CCCCACGT", "GGGGACGT", "TTTTACGT"]
    @test d.counts == [2, 2, 1, 1]
    @test [d.sequences[i] for i in d.map] == [r.sequence for r in reads]
    # Chunks of two reads: TTTT and AAAA are tabulated before CCCC and GGGG.
    d2 = dereplicate(reads; n=2)
    @test d2.sequences == ["AAAAACGT", "CCCCACGT", "TTTTACGT", "GGGGACGT"]
    @test [d2.sequences[i] for i in d2.map] == [r.sequence for r in reads]
    @test_throws ArgumentError dereplicate(reads; n=0)
end

@testset "Compact mean qualities" begin
    reads = [FastqRecord("a", "ACGTACGT", [30, 31, 32, 33, 34, 35, 36, 37]),
             FastqRecord("b", "ACGTACGT", [31, 31, 31, 31, 31, 31, 31, 40]),
             FastqRecord("c", "TTTTACGT", fill(20, 8))]
    d = dereplicate(reads)
    @test d.qualities isa DADA2.MeanQualities
    @test d.qualities[1] == ([30, 31, 32, 33, 34, 35, 36, 37] .+ [31, 31, 31, 31, 31, 31, 31, 40]) ./ 2
    @test d.qualities[2] == fill(20.0, 8)
    @test d.qualities == [[30.5, 31.0, 31.5, 32.0, 32.5, 33.0, 33.5, 38.5], fill(20.0, 8)]
    @test DADA2.rounded_qualities(d.qualities, 1) == UInt8[31, 31, 32, 32, 33, 33, 34, 39]
    @test DADA2.max_quality(d.qualities) == 38.5
end

@testset "Error learning with pooling identical to R (DD13)" begin
    samples = [uniques_derep(joinpath(PIPELINE, "pool$(s)_uniques.tsv")) for s in 1:4]
    matrix(path) = permutedims(hcat([parse.(Float64, r) for r in read_tsv(path)]...))
    for (mode, pool) in (("independent", false), ("pooled", true), ("pseudo", :pseudo)), threaded in (false, true)
        learned = learnErrors(samples; pool, threaded)
        # With pool=TRUE, R reports the pooled counts once per sample (each
        # sample's dada object carries them); the error model uses them once.
        copies = pool === true ? length(samples) : 1
        @test copies .* learned.transitions == round.(Int, matrix(joinpath(PIPELINE, "learn_$(mode)_trans.tsv")))
        @test learned.errors.rates ≈ matrix(joinpath(PIPELINE, "learn_$(mode)_err.tsv")) rtol=1e-12
        @test learned.iterations == parse(Int, readline(joinpath(PIPELINE, "learn_$(mode)_rounds.txt")))
    end
    root = joinpath(@__DIR__, "fixtures", "learn")
    files = [joinpath(root, "sample$i.fastq.gz") for i in 1:3]
    for mode in ("pooled", "pseudo")
        learned = learnErrors(files; pool=mode == "pooled" ? true : :pseudo)
        copies = mode == "pooled" ? length(files) : 1
        @test copies .* learned.transitions == round.(Int, matrix(joinpath(root, "$(mode)_trans.tsv")))
        @test learned.errors.rates ≈ matrix(joinpath(root, "$(mode)_err.tsv")) rtol=1e-12
    end
    @test_throws ArgumentError learnErrors(files; pool=true, qualityMode=:raw)
    @test_throws ArgumentError learnErrors(samples; pool=:unknown)
    @test learnErrors(files; pool=:pseudo, qualityMode=:raw).quality_mode == :raw
end

@testset "OMEGA bounds as R (DD12)" begin
    # Two 16-base sequences at Q40, 1,000 and 100 reads.
    parent = "ACGTACGTTGCATGCA"; variant = "ACGTACGTTGCATGCT"
    d = Derep([parent, variant], [1000, 100], [fill(40.0, 16), fill(40.0, 16)])
    @test length(dada(d).sequences) == 2
    @test dada(d; omegaA=0).sequences == [parent]   # no abundance births
    @test dada(d; omegaA=0).rejected == 100          # the variant fails the final test
    @test length(dada(d; omegaP=0).sequences) == 2
    @test dada(d; omegaA=0, omegaP=0).sequences == [parent]
    @test_throws ArgumentError dada(d; omegaA=1)
    @test_throws ArgumentError dada(d; omegaP=1)
    @test_throws ArgumentError dada(d; omegaA=-1e-3)
end
