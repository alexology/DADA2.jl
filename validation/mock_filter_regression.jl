using DADA2, Test

# Requires the public raw reads and R-filtered outputs from scripts/scientific.ps1.
# Use temporary outputs to preserve the frozen 0.2.0 scientific artifacts.
root = joinpath(@__DIR__, "generated", "public")
function records(path)
    reader = FastqReader(path)
    try
        collect(reader)
    finally
        close(reader)
    end
end
identical(a, b) = a.id == b.id && a.sequence == b.sequence && a.quality == b.quality

mktempdir() do dir
    f = joinpath(dir, "forward.fastq.gz")
    r = joinpath(dir, "reverse.fastq.gz")
    stats = filterAndTrim(joinpath(root, "Mock_S280_L001_R1_001.fastq"),
        joinpath(root, "Mock_S280_L001_R2_001.fastq"), f, r;
        optionsF=FilterOptions(truncLen=240,maxEE=2.0),
        optionsR=FilterOptions(truncLen=160,maxEE=2.0))
    jf = records(f); jr = records(r)
    rf = records(joinpath(root, "r_filteredF.fastq.gz"))
    rr = records(joinpath(root, "r_filteredR.fastq.gz"))
    @testset "Public mock filter parity with DADA2 R 1.36.0" begin
        @test (stats.reads_in, stats.reads_out) == (4779, 4314)
        @test length(jf) == length(jr) == length(rf) == length(rr) == 4314
        @test all(length(x.sequence) == 240 for x in jf)
        @test all(length(x.sequence) == 160 for x in jr)
        @test length(jf) == length(rf) && all(identical.(jf, rf))
        @test length(jr) == length(rr) && all(identical.(jr, rr))
    end
    path = joinpath(@__DIR__, "mock_filter_regression.tsv")
    open(path, "w") do io
        println(io, "julia_package_version\treads_in\tjulia_pairs\tr_pairs\tidentical_forward\tidentical_reverse")
        println(io, pkgversion(DADA2), '\t', stats.reads_in, '\t', length(jf), '\t', length(rf),
            '\t', count(identity, identical.(jf, rf)), '\t', count(identity, identical.(jr, rr)))
    end
end
