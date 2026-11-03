# DADA2.jl

Native Julia implementation of the [DADA2](https://github.com/benjjneb/dada2)
pipeline for amplicon sequence variants (ASVs): filtering, error learning,
denoising, pair merging, chimera removal and taxonomic assignment.

The denoising core, error learning, pair merging and chimera detection are
translations of the DADA2 1.36.0 R/C++ code. For the same inputs and options
they give the same results as R (see [Fidelity to DADA2 R](#fidelity-to-dada2-r)).
Taxonomic assignment is an independent implementation, not validated against
R. Version **0.3.0** plus unreleased changes ([CHANGELOG](CHANGELOG.md)).

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.
High performance is a goal, not a measured guarantee: see
[Performance](#performance) for what has actually been measured.

## Installation

The package is not registered in the Julia General registry. Install it from
GitHub (Julia ≥ 1.10; verified locally with Julia 1.12.2 on Windows):

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/DADA2.jl")
```

Start Julia with several threads (`julia --threads=4`, or `auto`) to use them.

## Paired-end pipeline

```julia
using DADA2

stats = filterAndTrim("R1.fastq.gz", "R2.fastq.gz", "F.fastq.gz", "R.fastq.gz";
    optionsF=FilterOptions(trimLeft=20, truncLen=240, maxEE=2.0),
    optionsR=FilterOptions(trimLeft=20, truncLen=160, maxEE=2.0))
errF = learnErrors("F.fastq.gz")
errR = learnErrors("R.fastq.gz")
dadaF = dada(derepFastq("F.fastq.gz"), errF)
dadaR = dada(derepFastq("R.fastq.gz"), errR)
merged = mergePairs(dadaF, dadaR)
table = removeBimeraDenovo(makeSequenceTable([merged]; samples=["sample1"]))
writeFasta("asvs.fasta", table.sequences)
```

The R equivalents are `filterAndTrim(..., trimLeft=c(20,20), truncLen=c(240,160),
maxEE=c(2,2))`, `learnErrors`, `derepFastq`, `dada`, `mergePairs`,
`makeSequenceTable` and `removeBimeraDenovo(method="consensus")`, with the same
defaults (including `truncQ=2`, `maxN=0`, `minLen=20` and phiX removal).
Paired filtering checks the read identifiers and keeps a pair only if both
reads pass. `mergePairs` needs the read maps kept by `derepFastq` (default).
Output files are never overwritten.

For single-end data, use `filterAndTrim(input, output; options=...)` and build
the table from the `dada` results. [examples/real_data.jl](examples/real_data.jl)
runs the whole pipeline on one or many samples from the command line
([instructions](docs/real-data-example.md)).

## Main functions

| Function | Behaviour |
| --- | --- |
| `filterAndTrim`, `FilterOptions`, `isPhiX` | R `filterAndTrim`: `trimLeft`, `truncLen`, `truncQ`, `maxN`, `maxEE`, `minLen`, `maxLen`, phiX removal (`rmPhix`) |
| `derepFastq`, `dereplicate` | R `derepFastq`: counts, mean qualities, read map; uniques in R's order |
| `learnErrors`, `errorDiagnostics` | R `learnErrors` with `loessErrfun`; training support by base and quality |
| `dada` | R `dada`: options of `setDadaOpt`, priors, `pool=true` and `pool=:pseudo`, `$clustering` diagnostics |
| `mergePairs` | R `mergePairs`: `minOverlap`, `maxMismatch`, `trimOverhang`, `justConcatenate`; all pairs in `pairs` |
| `makeSequenceTable` | Samples × ASVs table, columns by decreasing abundance |
| `removeBimeraDenovo`, `isBimeraDenovo`, `isBimeraDenovoTable`, `isBimera` | R chimera detection: `:consensus`, `:pooled`, `:per_sample`, `allowOneOff` |
| `trainTaxonomy`, `assignTaxonomy` | k-mer naive Bayes classifier with bootstrap; reusable index |
| `trainSpecies`, `addSpecies` | Exact matching of whole sequences or fragments; reusable index |
| `readFasta`, `writeFasta`, `FastqReader` | Plain or gzip FASTA/FASTQ (FASTX.jl) |

### Results

- `dada` returns a `DadaResult`: `sequences` and `counts` of the ASVs by
  decreasing abundance (R keeps the order of discovery), `assignments[i]`, the
  ASV of unique sequence `i` (zero if not corrected by the final test),
  `rejected` reads, `transitions`, `converged` and `clustering`, R's
  `$clustering` table in the order of `sequences`.
  `sum(result.counts) + result.rejected == sum(result.derep.counts)`.
- `learnErrors` returns an `ErrorLearningResult`; use `result.errors` (an
  `ErrorModel`, 16 transitions × qualities from Q0) and check
  `result.converged`. `dada(derep, learned)` accepts it directly.
- `mergePairs` returns a `MergeResult` with the merged `sequences` and
  `counts` (identical sequences summed), `accepted` and `rejected` reads, and
  `pairs`, R's data frame with `returnRejects=TRUE`.

### Options

`dada` accepts the options of R `setDadaOpt` under Julia names and with the
same defaults: `omegaA=1e-40`, `omegaP=1e-4`, `omegaC=1e-40`,
`detectSingletons`, `band=16`, `match=5`, `mismatch=-4`, `gap=-8`,
`homopolymerGap`, `useKmers`, `kdistCutoff=0.42`, `vectorizedAlignment`,
`gapless`, `greedy`, `maxClusters`, `minFold`, `minHamming`, `minAbundance`,
plus `priors`. With several samples, `pool=true` denoises the pooled uniques
once and `pool=:pseudo` repeats inference with the ASVs found in at least
`pseudoPrevalence=2` samples as priors.

`learnErrors` reads whole files until it has more than `nbases`
(100,000,000) bases and runs at most `maxIter=10` rounds, as R. For FASTQ
inputs, `qualityMode=:raw` counts transitions at the original read qualities
instead of the mean qualities of each unique sequence: this is a Julia
extension, not validated (see [Fidelity](#fidelity-to-dada2-r)).

`phredErrors()` gives a theoretical error model for tests; it is not
calibrated to any instrument. A model estimated elsewhere can be imported as
`ErrorModel(matrix)` (rows A→A, A→C, …, T→T, each block of four summing to
one; columns from Q0).

Sequences can be given as strings or BioSequences `LongDNA` values.
Dereplicated sequences must contain only A/C/G/T and be longer than 5 bases,
as in R. Multi-line FASTQ and bzip2 are not supported.

### Taxonomy

```julia
taxonomy_model = trainTaxonomy("training.fa.gz"; k=8)
taxonomy = assignTaxonomy(table.sequences, taxonomy_model; seed=42)
species_model = trainSpecies("species.fa.gz")
species = addSpecies(table.sequences, species_model)
```

Training headers must contain contiguous semicolon-separated ranks; empty
internal ranks are rejected. References may contain IUPAC codes (k-mers
spanning them are skipped); species matches are literal.

## Fidelity to DADA2 R

With the same inputs and options, the following give R's results, checked
against R dada2 1.36.0 (R 4.5.0):

- `filterAndTrim`, `derepFastq`, `dada` (all options, pooling, diagnostics),
  `learnErrors`, `mergePairs` and the chimera functions, on fixtures generated
  by R scripts in `validation/` and kept in `test/fixtures/`;
- the whole pipeline on a real paired-end COI sample (292,050 read pairs,
  Element AVITI): identical filtered reads, error models within 1.3e-14,
  identical ASVs, merged pairs and final table (303 ASVs).

Floating-point results agree to rounding (about 1e-14 for error models): R's
`loess` uses LINPACK and `ppois` R's own Mathlib, Julia LAPACK and
Distributions.jl. Deliberate differences and what is missing are listed in
[docs/r-parity.md](docs/r-parity.md), with the commands to repeat the checks.
In short:

- sequences above about 6,500 bases are aligned with 32-bit scores, where
  R's 16-bit scores overflow;
- an error model that does not cover the qualities is rejected rather than
  extended;
- not ported: primer removal (`removePrimers`), low-complexity filtering,
  PacBio mode, `collapseNoMismatch`, plots and some other helpers;
  taxonomic assignment is not R's RDP implementation.

A [preliminary scientific validation](docs/scientific-validation.md) of
version 0.2.0 compared an earlier, independent algorithm with R; it does not
describe the current code.

## Performance

On the real paired-end sample above (Windows, Intel Core i7-9750H, four
threads; one run each, Julia compilation excluded):

| Step | R | Julia |
| --- | ---: | ---: |
| Filtering | 44 s | 19 s |
| Error learning | 212 s | 47 s |
| Dereplication | 5 s | 2 s |
| Denoising | 38 s | 21 s |
| Merging | 2 s | 2 s |
| Chimeras | 10 s | 4 s |
| **Total** | **312 s** | **95 s** |

The main gains come from keeping alignments across error-learning rounds
(R recomputes them) and from a vectorized k-mer screen. This is one sample on
one machine. The first call in a new Julia session also compiles the code
(tens of seconds), and on small samples with one thread R can be faster (0.51
s against 0.71 s for 4,152 unique sequences). Repeat the comparison on your
data with [benchmark/pipeline_timing.jl](benchmark/pipeline_timing.jl) and
[benchmark/pipeline_timing.R](benchmark/pipeline_timing.R).

Several threads are used for the comparisons inside `dada`, for samples in
parallel and for the preparation of the uniques. The alignment cache speeds up
error learning; `cacheAlignments=false` saves memory. The k-mer counts take
1 KB per unique sequence. On several nodes, learn the error model once and
distribute the per-sample inference through the job scheduler.

## Development

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 test/runtests.jl
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the R comparisons that changes to
scientific behaviour must pass.

Licensed under LGPL-3.0-only; references and provenance in [NOTICE.md](NOTICE.md).
