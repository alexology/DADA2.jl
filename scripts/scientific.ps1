param([string]$RscriptPath = 'Rscript', [switch]$SkipDownload)
$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
Push-Location $projectRoot
try {
    $publicRoot = Join-Path $projectRoot 'validation/generated/public'
    $null = New-Item -ItemType Directory -Path $publicRoot -Force
    $zipPath = Join-Path $publicRoot 'miseqsopdata.zip'
    if (!$SkipDownload -and !(Test-Path -LiteralPath $zipPath)) {
        Invoke-WebRequest -Uri 'https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip' -OutFile $zipPath
    }
    $expected = '728f7544824f8ca0223618b18bb60fce6a5dc20ae3d476267dbd5cf5594a35d4'
    if ((Get-FileHash -LiteralPath $zipPath -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected) {
        throw 'Public dataset checksum mismatch'
    }
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip = [IO.Compression.ZipFile]::OpenRead($zipPath)
    try {
        foreach ($entry in $zip.Entries) {
            if ($entry.Name -in @('Mock_S280_L001_R1_001.fastq','Mock_S280_L001_R2_001.fastq','HMP_MOCK.v35.fasta')) {
                $destination = Join-Path $publicRoot $entry.Name
                if (!(Test-Path -LiteralPath $destination)) {
                    [IO.Compression.ZipFileExtensions]::ExtractToFile($entry,$destination,$false)
                }
            }
        }
    } finally { $zip.Dispose() }
    & julia --project=. --threads=4 validation/scientific_generate.jl
    if ($LASTEXITCODE) { throw 'Simulation generation failed' }
    & $RscriptPath validation/scientific_reference.R
    if ($LASTEXITCODE) { throw 'R simulation reference failed' }
    & julia --project=. --threads=4 validation/scientific_evaluate.jl
    if ($LASTEXITCODE) { throw 'Simulation evaluation failed' }
    & $RscriptPath validation/mock_reference.R
    if ($LASTEXITCODE) { throw 'R mock reference failed' }
    & julia --project=. --threads=4 validation/mock_evaluate.jl
    if ($LASTEXITCODE) { throw 'Experimental mock evaluation failed' }
    & julia --project=. validation/mock_filter_audit.jl
    if ($LASTEXITCODE) { throw 'Filter audit failed' }
    & $RscriptPath validation/scientific_summary.R
    if ($LASTEXITCODE) { throw 'Scientific summary failed' }
    & (Join-Path $PSScriptRoot 'scientific-metadata.ps1')
} finally { Pop-Location }
