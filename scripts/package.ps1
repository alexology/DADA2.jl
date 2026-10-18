param([string]$OutputDirectory = "")
$ErrorActionPreference = 'Stop'
$projectRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
if (!$OutputDirectory) { $OutputDirectory = Join-Path $projectRoot 'dist' }
$OutputDirectory = [IO.Path]::GetFullPath($OutputDirectory)
$versionLine = Get-Content -LiteralPath (Join-Path $projectRoot 'Project.toml') | Where-Object { $_ -match '^version\s*=' }
$packageVersion = [regex]::Match($versionLine, '"([0-9]+\.[0-9]+\.[0-9]+)"').Groups[1].Value
if (!$packageVersion) { throw 'Missing semantic package version' }
$null = New-Item -ItemType Directory -Path $OutputDirectory -Force
$archivePath = Join-Path $OutputDirectory "DADA2.jl-$packageVersion.zip"
if (Test-Path -LiteralPath $archivePath) { throw "Archive already exists: $archivePath" }
$temporaryPath = "$archivePath.tmp"
if (Test-Path -LiteralPath $temporaryPath) { throw "Temporary archive already exists: $temporaryPath" }
$directories = @('src', 'test', 'examples', 'docs', 'benchmark', 'validation', 'scripts', '.github')
$files = @('Project.toml', 'README.md', 'LICENSE', 'NOTICE.md', 'CHANGELOG.md', 'CONTRIBUTING.md', '.gitignore')
$included = @()
foreach ($name in $files) { $included += Get-Item -LiteralPath (Join-Path $projectRoot $name) }
foreach ($name in $directories) {
    $included += Get-ChildItem -LiteralPath (Join-Path $projectRoot $name) -Recurse -File -Force | Where-Object {
        $_.FullName -notlike '*\validation\generated\*' -and
        $_.FullName -notlike '*/validation/generated/*' -and $_.Name -ne 'results.tsv'
    }
}
Add-Type -AssemblyName System.IO.Compression
Add-Type -AssemblyName System.IO.Compression.FileSystem
$archive = [IO.Compression.ZipFile]::Open($temporaryPath, [IO.Compression.ZipArchiveMode]::Create)
try {
    foreach ($file in ($included | Sort-Object FullName)) {
        $relative = $file.FullName.Substring($projectRoot.Length + 1).Replace('\','/')
        $entry = $archive.CreateEntry("DADA2.jl/$relative", [IO.Compression.CompressionLevel]::Optimal)
        $entry.LastWriteTime = [DateTimeOffset]::new(2020,1,1,0,0,0,[TimeSpan]::Zero)
        $source = [IO.File]::OpenRead($file.FullName)
        try {
            $destination = $entry.Open()
            try { $source.CopyTo($destination) } finally { $destination.Dispose() }
        } finally { $source.Dispose() }
    }
} catch {
    $archive.Dispose()
    Remove-Item -LiteralPath $temporaryPath
    throw
} finally { $archive.Dispose() }
[IO.File]::Move($temporaryPath,$archivePath)
$hash = (Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash.ToLowerInvariant()
[IO.File]::WriteAllText("$archivePath.sha256", "$hash  $([IO.Path]::GetFileName($archivePath))`n")
Write-Output "$archivePath`nSHA256: $hash`nFiles: $($included.Count)"
