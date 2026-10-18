$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$sourceFiles = Get-ChildItem -LiteralPath (Join-Path $projectRoot 'src') -File | Sort-Object Name
$sourceHashes = @()
foreach ($sourceFile in $sourceFiles) {
    $sourceHashes += @{file="src/$($sourceFile.Name)"; sha256=(Get-FileHash -LiteralPath $sourceFile.FullName -Algorithm SHA256).Hash.ToLowerInvariant()}
}
Write-Output 'Source hashes recorded'
$sampleRoot = Join-Path $projectRoot 'validation/generated/science'
$samples = @()
$caseDirectories = Get-ChildItem -LiteralPath $sampleRoot -Directory | Sort-Object Name
foreach ($caseDirectory in $caseDirectories) {
    $hashes = @{}
    foreach ($name in @('reads.fastq','truth.tsv','oracle.tsv')) {
        $hashes[$name] = (Get-FileHash -LiteralPath (Join-Path $caseDirectory.FullName $name) -Algorithm SHA256).Hash.ToLowerInvariant()
    }
    $samples += @{case=$caseDirectory.Name; hashes=$hashes}
}
Write-Output 'Simulated input hashes recorded'
$mockRoot = Join-Path $projectRoot 'validation/generated/public'
$publicHashes = @{}
foreach ($name in @('miseqsopdata.zip','HMP_MOCK.v35.fasta','Mock_S280_L001_R1_001.fastq','Mock_S280_L001_R2_001.fastq')) {
    $publicHashes[$name] = (Get-FileHash -LiteralPath (Join-Path $mockRoot $name) -Algorithm SHA256).Hash.ToLowerInvariant()
}
$juliaVersion = & julia --version
Write-Output 'Public data hashes and Julia version recorded'
if ($LASTEXITCODE) { throw 'Could not determine Julia version' }
$metadata = @{
    julia=$juliaVersion.ToString(); r_versions=[IO.File]::ReadAllLines((Join-Path $sampleRoot 'r_versions.txt'));
    denoiser_version='0.2.0'; source_hashes=$sourceHashes; simulated_inputs=$samples;
    public_input_hashes=$publicHashes; seeds=@(11,23,37,51,71);
    independent_communities=5; conditions=30; measured_method_sample_pairs=150
}
[IO.File]::WriteAllText((Join-Path $projectRoot 'validation/scientific_metadata.json'),
    ($metadata | ConvertTo-Json -Depth 8),[Text.Encoding]::UTF8)
Write-Output 'Scientific source and input hashes recorded'
