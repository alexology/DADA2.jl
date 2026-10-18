function fixture_matrix(path, type)
    rows = [parse.(type,split(line,'\t')) for line in readlines(path)]
    permutedims(hcat(rows...))
end

@testset "Error learning: smoothing and training support" begin
    root = joinpath(@__DIR__,"fixtures")
    trans = fixture_matrix(joinpath(root,"error_transitions.tsv"),Int)
    reference = fixture_matrix(joinpath(root,"error_loess_direct.tsv"),Float64)
    initial = phredErrors(maxQ=50)
    fitted = DADA2.fitted_errors(trans,initial)
    @test fitted.rates ≈ reference rtol=1e-10 atol=1e-13
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
    a = fixture([repeat("A",40),"C"*repeat("A",39)],[1000,2];q=30)
    b = fixture([repeat("C",40),"A"*repeat("C",39)],[2000,2];q=40)
    learned = learnErrors([a,b];maxIter=1,threaded=false)
    results = dada([a,b],phredErrors();threaded=false,omegaC=0)
    @test learned.transitions == results[1].transitions + results[2].transitions
    @test learned.errors.rates == DADA2.fitted_errors(learned.transitions,phredErrors()).rates
    threaded = learnErrors([a,b];maxIter=1,threaded=true)
    @test threaded.transitions == learned.transitions
    @test threaded.errors.rates == learned.errors.rates
    low_memory = learnErrors([a,b];maxIter=1,threaded=true,cacheAlignments=false)
    @test low_memory.transitions == learned.transitions
    @test low_memory.errors.rates == learned.errors.rates
    diag = errorDiagnostics(learned)
    @test sum(diag.exposures) == sum(learned.transitions)
    @test sum(diag.observed_substitutions) == sum(learned.transitions[[2,3,4,5,7,8,9,10,12,13,14,15],:])
    @test diag.supported == (diag.exposures .> 0)
    @test diag.unobserved_sources == ['G','T']
    @test diag.expected_substitutions[.!diag.supported] == zeros(count(!,diag.supported))
    @test diag.inference_converged == all(r.converged for r in results)
    @test_throws ArgumentError learnErrors(a;smoother=:unknown)
    @test_throws ArgumentError learnErrors(a;tolerance=Inf)
    @test_throws ArgumentError learnErrors(a;qualityMode=:raw)
    @test learnErrors(a;smoother=:linear,maxIter=1).errors.rates ==
        DADA2.linear_fitted_errors(dada(a).transitions,phredErrors()).rates
    mktempdir() do dir
        path = joinpath(dir,"mixed.fastq")
        # Same sequence, different qualities: mean counts collapse to Q30;
        # raw counts must retain ten base observations each at Q20 and Q40.
        write(path,"@a\nAAAAAAAAAA\n+\n5555555555\n@b\nAAAAAAAAAA\n+\nIIIIIIIIII\n")
        raw = learnErrors(path;maxIter=1)
        averaged = learnErrors(path;qualityMode=:mean,maxIter=1)
        @test raw.quality_mode == :raw
        @test averaged.quality_mode == :mean
        @test raw.transitions[1,21] == raw.transitions[1,41] == 10
        @test sum(raw.transitions) == 20
        @test raw.training_reads == 2
        @test raw.training_bases == 20
        @test averaged.transitions[1,31] == 20
        @test averaged.errors.rates == learnErrors(derepFastq(path);maxIter=1).errors.rates
        @test errorDiagnostics(raw).supported[1,21] && errorDiagnostics(raw).supported[1,41]
        @test !errorDiagnostics(raw).supported[1,31]
        @test errorDiagnostics(raw).training_reads == 2
        @test errorDiagnostics(raw).training_bases == 20
        @test_throws ArgumentError learnErrors(path;qualityMode=:unknown)
        @test_throws ArgumentError learnErrors(path;nbases=0)
        @test_throws ArgumentError learnErrors(String[])
        # Mean Q30 is covered, but the original Q40 is outside this model.
        @test_throws ArgumentError learnErrors(path;initial=phredErrors(maxQ=35),maxIter=1)
        pooled = learnErrors([path,path];maxIter=1)
        @test pooled.transitions == 2raw.transitions
        # The budget is divided deterministically between files. Whole reads
        # may cross the exact base limit, but both equal-size inputs contribute.
        budgeted = learnErrors([path,path];nbases=20,maxIter=1,threaded=false)
        budgeted_threaded = learnErrors([path,path];nbases=20,maxIter=1,threaded=true)
        @test budgeted.training_reads == 2
        @test budgeted.training_bases == 20
        @test budgeted.transitions[1,21] == 20
        @test sum(budgeted.transitions) == 20
        @test budgeted_threaded.transitions == budgeted.transitions
        @test budgeted_threaded.errors.rates == budgeted.errors.rates
        crossed = learnErrors(path;nbases=11,maxIter=1)
        @test crossed.training_bases == 20
        @test crossed.training_reads == 2
        noisy = joinpath(dir,"substitutions.fastq")
        write(noisy,join(["@parent$i\nACGTACGTACGT\n+\nIIIIIIIIIIII\n" for i in 1:100]) *
                    "@error1\nACGTTCGTACGT\n+\n555555555555\n" *
                    "@error2\nACGTTCGTACGT\n+\nIIIIIIIIIIII\n")
        substitutions = learnErrors(noisy;maxIter=1,useKmers=false)
        @test substitutions.transitions[4,21] == substitutions.transitions[4,41] == 1
        @test sum(substitutions.transitions) == 1224
        @test substitutions.transitions[1,41] == 302
    end
end
