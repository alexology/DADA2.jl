# Contributing

This package is an experimental Julia implementation of the central DADA2
pipeline. Read the documented differences from R before interpreting results.

Use a branch for each change. Run from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 -e 'using Pkg; Pkg.test()'
julia --project=. --threads=4 benchmark/benchmarks.jl
```

Changes to scientific behavior should include an independent reference check
or a regression fixture. The R comparison scripts in `validation/` compare
the same error model, keeping inference separate from error learning.
Do not claim general equivalence from a single fixture.

Report the Julia version, operating system, parameters and smallest input
that reproduces a bug. Include the complete exception and stack trace.
For performance changes, compare identical inputs after compilation and
report allocations and thread count as well as runtime.

## Source release

On Windows use `powershell -NoProfile -ExecutionPolicy Bypass -File scripts/package.ps1`
(the policy option only applies to that process); on Linux/macOS with
PowerShell installed use `pwsh -File scripts/package.ps1`. The script produces
`dist/DADA2.jl-<version>.zip` and its SHA-256 checksum and refuses to replace
an existing archive. Archives exclude the local Manifest, upstream checkout,
generated validation results, Git metadata and previous archives.

The GitHub **Source archive** workflow builds and checks an archive using
`workflow_dispatch`. Its artifact can be downloaded and attached to a release.
Version changes belong in Project.toml and CHANGELOG.md.
