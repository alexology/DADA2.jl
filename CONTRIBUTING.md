# Contributing

DADA2.jl translates DADA2 1.36.0: for the same inputs and options it must give
R's results. Read [docs/r-parity.md](docs/r-parity.md) before changing
filtering, dereplication, denoising, error learning, merging or chimera
detection.

Use a branch for each change. Run from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 test/runtests.jl
```

`Pkg.test()` runs the same suite with bounds checking (slower).

## Changes to scientific behaviour

- Keep upstream's iteration order, tie-breaking and order of floating-point
  operations where they can change a result; the tests compare exact ASVs,
  maps and transition counts with R.
- A deliberate difference from R must be documented in `docs/r-parity.md`
  and the CHANGELOG, with a test showing it.
- To extend the R comparisons, change the generators in `validation/`
  (`core_reference.R`, `learn_reference.R`, `pipeline_reference.R`,
  `phix_reference.R`) and rerun them with R and dada2 installed: they rewrite
  the fixtures in `test/fixtures/`. Do not edit fixtures by hand.
- Do not claim general equivalence from one dataset; for real data, compare
  the whole pipeline with `benchmark/pipeline_timing.{R,jl}` and
  `benchmark/pipeline_compare.jl`. Keep private data out of the repository.

## Performance changes

Results must not change: run the test suite and, for changes to the
denoising core, compare a real or simulated sample with R as above. Compare
identical inputs after compilation and report thread count, run time and
allocations; `benchmark/benchmarks.jl` gives quick micro-benchmarks.

## Bug reports

Report the Julia version, operating system, parameters and the smallest input
that reproduces the problem, with the complete exception and stack trace.

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
