param([string]$RscriptPath = 'Rscript')
$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
Push-Location $projectRoot
try {
    foreach ($name in @('Mock_S280_L001_R1_001.fastq','Mock_S280_L001_R2_001.fastq','mock_errF.tsv','mock_errR.tsv')) {
        if (!(Test-Path -LiteralPath (Join-Path 'validation/generated/public' $name))) { throw "Missing public input: $name" }
    }
    # Run each process sequentially: competing benchmarks distort measurements.
    foreach ($threadCount in @(1,4)) {
        & $RscriptPath --vanilla benchmark/mock_timing.R $threadCount "benchmark/mock_timing_R_$threadCount.tsv"
        if ($LASTEXITCODE) { throw "R benchmark failed ($threadCount threads)" }
        & julia --startup-file=no --project=. "--threads=$threadCount" benchmark/mock_timing.jl "benchmark/mock_timing_Julia_$threadCount.tsv"
        if ($LASTEXITCODE) { throw "Julia benchmark failed ($threadCount threads)" }
    }
} finally { Pop-Location }
