$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$processor = Get-CimInstance Win32_Processor | Select-Object -First 1
$operatingSystem = Get-CimInstance Win32_OperatingSystem
$relativePaths = @('benchmark/mock_timing.jl','benchmark/mock_timing.R','benchmark/mock_timing_summary.R',
    'scripts/mock-timing.ps1','src/inference.jl','src/alignment.jl','src/fastq.jl','src/community.jl',
    'Project.toml','validation/generated/public/Mock_S280_L001_R1_001.fastq',
    'validation/generated/public/Mock_S280_L001_R2_001.fastq','validation/generated/public/mock_errF.tsv',
    'validation/generated/public/mock_errR.tsv','validation/generated/public/r_mock.tsv')
$hashes = foreach ($relative in $relativePaths) {
    [ordered]@{path=$relative; sha256=(Get-FileHash -LiteralPath (Join-Path $projectRoot $relative) -Algorithm SHA256).Hash.ToLowerInvariant()}
}
$metadata = [ordered]@{
    cpu=$processor.Name.ToString(); physical_cores=$processor.NumberOfCores; logical_processors=$processor.NumberOfLogicalProcessors;
    operating_system=$operatingSystem.Caption.ToString(); os_version=$operatingSystem.Version.ToString();
    julia_version=(& julia --version).ToString();
    r_versions=[IO.File]::ReadAllLines((Join-Path $projectRoot 'benchmark/mock_timing_R_versions.txt'));
    thread_counts=@(1,4); warm_repetitions=5; first_calls=10; total_measurements=60;
    blas_thread_setting=1; execution_order=@('R_1','Julia_1','R_4','Julia_4');
    timing='Wall time; GC outside timed pipeline. First script timing excludes interpreter startup and package installation/precompilation.';
    hashes=@($hashes)
}
[IO.File]::WriteAllText((Join-Path $projectRoot 'benchmark/mock_timing_metadata.json'),(ConvertTo-Json -InputObject $metadata -Depth 5))
