function fixture_matrix(path, type)
    rows = [parse.(type,split(line,'\t')) for line in readlines(path)]
    permutedims(hcat(rows...))
end

@testset "Error learning: smoothing and training support" begin
    root = joinpath(@__DIR__,"fixtures")
    trans = fixture_matrix(joinpath(root,"error_transitions.tsv"),Int)
    reference = fixture_matrix(joinpath(root,"error_loess_direct.tsv"),Float64)
    initial = phredErrors(maxQ=50)
    fitted = DADA2.fitted_errors(trans,initial;smoother=:direct)
    @test fitted.rates ≈ reference rtol=1e-10 atol=1e-13
    # Default: DADA2 loessErrfun on R's interpolated surface (fixture: 15 digits).
    interpolated = fixture_matrix(joinpath(root,"error_loess_r.tsv"),Float64)
    @test DADA2.fitted_errors(trans,initial).rates ≈ interpolated rtol=1e-12
    @test DADA2.r_loess_errors(trans,nothing).rates ≈ interpolated rtol=1e-12
    @test_throws ArgumentError DADA2.r_loess_errors(zeros(Int,16,51),nothing)
    @test fitted.rates[:,1:10] == repeat(fitted.rates[:,11:11],1,10)
    @test fitted.rates[:,42:51] == repeat(fitted.rates[:,41:41],1,10)
    @test all(isfinite,fitted.rates)
    @test DADA2.fitted_errors(zeros(Int,16,51),initial).rates == initial.rates
    @test_throws ArgumentError DADA2.fitted_errors(trans,initial;span=NaN)
    @test_throws ArgumentError DADA2.fitted_errors(trans,initial;smoother=:unknown)
    sparse = zeros(Int,16,51)
    sparse[1,31] = 10000; sparse[2,31] = 10
    one = DADA2.fitted_errors(sparse,initial)
    @test one.rates[2,:] ≈ fill(11/10010,51)
    @test one.rates[5:16,:] == initial.rates[5:16,:]
    sparse[1,41] = 20000; sparse[2,41] = 5
    @test all(isfinite,DADA2.fitted_errors(sparse,initial).rates)
    @test DADA2.fitted_errors(sparse,initial).rates[2,36] ≈ sqrt((11/10010)*(6/20005))
    @test all(isfinite,DADA2.fitted_errors(sparse,initial;span=0.01).rates)
    # Learning follows R learnErrors (parity in test_core.jl); here the
    # counting modes and options, on simulated reads without indels.
    files = [joinpath(root,"learn","sample$i.fastq.gz") for i in 1:3]
    learned = learnErrors(files;threaded=false)
    @test learned.quality_mode == :mean
    @test (learned.training_reads,learned.training_bases) == (1224,122400)
    threaded = learnErrors(files;threaded=true)
    low_memory = learnErrors(files;threaded=true,cacheAlignments=false)
    @test threaded.transitions == learned.transitions == low_memory.transitions
    @test threaded.errors.rates == learned.errors.rates == low_memory.errors.rates
    # Mean mode: every base counted once, at the rounded mean quality of its unique.
    dereps = derepFastq.(files)
    mean_hist = zeros(Int,size(learned.transitions,2))
    for d in dereps, (q,c) in zip(d.qualities,d.counts), x in q
        mean_hist[round(Int,x,RoundNearestTiesAway)+1] += c
    end
    @test vec(sum(learned.transitions;dims=1)) == mean_hist
    # Raw mode: every base counted at its original quality.
    raw = learnErrors(files;qualityMode=:raw)
    @test raw.quality_mode == :raw
    raw_hist = zeros(Int,size(raw.transitions,2))
    for path in files
        reader = FastqReader(path)
        foreach(read -> foreach(q -> raw_hist[q+1] += 1, read.quality), reader)
        close(reader)
    end
    @test vec(sum(raw.transitions;dims=1)) == raw_hist
    @test learnErrors(files;qualityMode=:raw,threaded=false).errors.rates == raw.errors.rates
    @test raw.errors.rates != learned.errors.rates
    @test all(q -> all(b -> sum(raw.errors.rates[4b+1:4b+4,q]) ≈ 1,0:3),axes(raw.errors.rates,2))
    # The same R sampling: whole files until the bases exceed nbases.
    budget = learnErrors(files;qualityMode=:raw,nbases=sum(length(r.quality) for r in FastqReader(files[1]))+1)
    @test budget.training_reads == sum(sum(d.counts) for d in dereps[1:2])
    diag = errorDiagnostics(learned)
    @test sum(diag.exposures) == sum(learned.transitions)
    @test sum(diag.observed_substitutions) == sum(learned.transitions[[2,3,4,5,7,8,9,10,12,13,14,15],:])
    @test diag.supported == (diag.exposures .> 0)
    @test isempty(diag.unobserved_sources)
    @test diag.expected_substitutions[.!diag.supported] == zeros(count(!,diag.supported))
    @test (diag.training_reads,diag.training_bases) == (1224,122400)
    for option in ((initial=phredErrors(),),(tolerance=1e-3,),(smoother=:linear,),(upstream=false,))
        @test_throws ArgumentError learnErrors(files;option...)
    end
    @test_throws ArgumentError learnErrors(dereps;qualityMode=:raw)
    @test_throws ArgumentError learnErrors(files;qualityMode=:unknown)
    @test_throws ArgumentError learnErrors(files;nbases=0)
    @test_throws ArgumentError learnErrors(String[])
    mktempdir() do dir
        # As in R, a base never observed leaves the rates undefined.
        path = joinpath(dir,"only_a.fastq")
        write(path,"@a
AAAAAAAAAA
+
5555555555
@b
AAAAAAAAAA
+
IIIIIIIIII
")
        @test_throws ArgumentError learnErrors(path)
    end
end
