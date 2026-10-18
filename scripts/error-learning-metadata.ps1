$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$paths = @('Project.toml','src/DADA2.jl','src/inference.jl','src/fastq.jl',
    'src/types.jl','src/alignment.jl','src/community.jl','src/taxonomy.jl',
    'validation/error_learning_evaluate.jl','validation/error_learning_mock.jl',
    'validation/error_learning_summary.R','validation/error_smoothing_fixture.jl',
    'validation/error_smoothing_reference.R','test/error_learning.jl',
    'scripts/error-learning-metadata.ps1','examples/pipeline.jl')
$sourceHashes = foreach ($relative in $paths) {
    $path = Join-Path $projectRoot $relative
    [ordered]@{path=$relative; sha256=(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()}
}
$manifest = Import-Csv -LiteralPath (Join-Path $projectRoot 'validation/generated/science/manifest.tsv') -Delimiter "`t"
$inputHashes = foreach ($case in $manifest) {
    foreach ($name in @('reads.fastq','truth.tsv','oracle.tsv','r_errors.tsv','r_learned.tsv')) {
        $relative = "validation/generated/science/$($case.case)/$name"
        [ordered]@{path=$relative; sha256=(Get-FileHash -LiteralPath (Join-Path $projectRoot $relative) -Algorithm SHA256).Hash.ToLowerInvariant()}
    }
}
$mockHashes = foreach ($name in @('r_filteredF.fastq.gz','r_filteredR.fastq.gz','HMP_MOCK.v35.fasta')) {
    $relative = "validation/generated/public/$name"
    [ordered]@{path=$relative; sha256=(Get-FileHash -LiteralPath (Join-Path $projectRoot $relative) -Algorithm SHA256).Hash.ToLowerInvariant()}
}
$packageVersion = [regex]::Match([IO.File]::ReadAllText((Join-Path $projectRoot 'Project.toml')),'(?m)^version\s*=\s*"([^"]+)"').Groups[1].Value
$metadata = [ordered]@{
    package_version=$packageVersion; julia_version=(& julia --version).ToString();
    r_versions=[IO.File]::ReadAllLines((Join-Path $projectRoot 'test/fixtures/error_loess_versions.txt'));
    seeds=@(11,23,37,51,71); cases=30; methods=7; measurements=210;
    source_hashes=@($sourceHashes); simulation_input_hashes=@($inputHashes); mock_input_hashes=@($mockHashes)
}
$json = ConvertTo-Json -InputObject $metadata -Depth 6
[IO.File]::WriteAllText((Join-Path $projectRoot 'validation/error_learning_metadata.json'),$json)
