# DADA2.jl

Native Julia implementation of the core [DADA2](https://github.com/benjjneb/dada2)
pipeline, with streaming FASTQ reading, statistical ASV inference and optimised
sequence comparisons.
Version **0.3.0**: it is not yet a complete port, nor a scientifically
validated replacement for all the functionality of the R package.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

A [preliminary scientific validation](docs/scientific-validation.md) is now
available: a real mock community and 30 simulated conditions with known
ground truth, compared with DADA2 R. The mock yields 20 ASVs, all matching the
references; the simulations reveal missed rare variants and compatibility
differences. Protocol, individual results and plots are included in the project.
The scientific report refers to 0.2.0; 0.2.1 fixes the `truncLen` filter, with
[regression tests and a comparison with R](docs/release-0.2.1.md).
0.3.0 adds quadratic LOESS, counts at the original FASTQ qualities and
learning diagnostics: [method and results](docs/error-learning.md).
A [timing comparison with DADA2 R on the mock](docs/mock-timing.md) is also
available: on the complete pipeline, R is faster in this local test.

## Install from GitHub

This package is not registered in the Julia General registry. Install it
from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alexology/DADA2.jl")
```

## Installation and usage

Requires Julia ≥ 1.10. Verified locally with Julia 1.12.2 on Windows.
From this directory:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 test/runtests.jl
julia --project=. --threads=4 examples/pipeline.jl reads.fastq results
```

To use the package from another project:

```julia
using Pkg
Pkg.develop(path="D:/r_packages/dada2_julia")
using DADA2, BioSequences

alignment = alignSequences(dna"ACGTACGT", dna"ACGTTCGT")
```

## Single-end pipeline

To try real single-end or paired-end data, including multiple samples, use
[examples/real_data.jl](examples/real_data.jl): [instructions and parameters](docs/real-data-example.md).

```julia
using DADA2

stats = filterAndTrim("sample.fastq.gz", "filtered.fastq.gz";
    options=FilterOptions(truncLen=240, maxEE=2.0, minLen=20))
reads = derepFastq("filtered.fastq.gz"; keepmap=false)
learned = learnErrors("filtered.fastq.gz"; nbases=100_000_000)
result = dada(reads, learned)
table = removeBimeraDenovo(makeSequenceTable([result]; samples=["sample"]))
writeFasta("asvs.fasta", table.sequences)

# Check convergence before interpreting the results.
@show learned.converged result.converged result.rejected
```

`phredErrors()` provides a theoretical starting model, not a model calibrated
to your instrument. A previously estimated model can be imported as
`ErrorModel(matrix)`: 16 rows A→A, A→C, …, T→T and one column per Q starting
from Q0. Each block of four transitions must sum to one.
Mean qualities are rounded half up, as in the upstream kernels.

## Paired-end pipeline

```julia
filterAndTrim("F.fastq.gz", "R.fastq.gz", "filteredF.fastq.gz", "filteredR.fastq.gz";
    optionsF=FilterOptions(truncLen=240, maxEE=2.0),
    optionsR=FilterOptions(truncLen=160, maxEE=2.0))
f = derepFastq("filteredF.fastq.gz")
r = derepFastq("filteredR.fastq.gz")
df = dada(f, learnErrors("filteredF.fastq.gz"))
dr = dada(r, learnErrors("filteredR.fastq.gz"))
merged = mergePairs(df, dr; minOverlap=12, maxMismatch=0)
table = removeBimeraDenovo(makeSequenceTable([merged]))
```

Filtering checks the number and identifiers of the paired reads and keeps the
pairs in which both mates are valid. `mergePairs` requires the read map: keep
`keepmap=true`, the default. If you dereplicate files that were filtered
separately, the caller is responsible for the order and identity of the pairs.
Existing output files are not overwritten.

## Features

| API | Implemented behaviour |
| --- | --- |
| `FastqReader` | Four-line FASTQ iterator, plain and gzip; close with `close` |
| `filterAndTrim` | Left trimming, truncation by length and quality, maxN/maxEE, length limits; synchronised paired-end |
| `derepFastq`, `dereplicate` | Counts, mean per-position qualities and optional read map |
| `dada` | Iterative partitioning, conditional Poisson, Bonferroni, priors, optional singletons, final OMEGA_C correction |
| `learnErrors`, `errorDiagnostics` | Weighted quadratic LOESS, raw counts or mean qualities, base budget, model shared across samples and estimate support |
| `mergePairs` | Suffix/prefix overlap without indels; accepted and rejected counts |
| `makeSequenceTable` | Samples × ASVs matrix with deterministic order |
| `isBimeraDenovo`, `removeBimeraDenovo` | Exact two-parent chimeras of equal length, pooled or consensus check |
| `trainTaxonomy`, `assignTaxonomy` | k-mer naive Bayes, reusable sparse index, IUPAC references, reproducible bootstrap |
| `trainSpecies`, `addSpecies` | Exact matching including fragments, reusable seed index, reverse complement and optional multi-species ambiguity |
| `readFasta`, `writeFasta` | Plain/gzip FASTA with IUPAC bases, optional BioSequences and empty-record check |

`DadaResult.assignments[i]` identifies the ASV of unique sequence `i`; zero
marks a read not corrected by the final check. Therefore
`sum(result.counts) + result.rejected == sum(result.derep.counts)`.
`omegaC=0` corrects all assigned sequences. `omegaA`, `omegaP` and `omegaC`
are the statistical thresholds for inference and correction.

For multiple samples, `dada(vector_of_dereps, errors)` parallelises over
samples and `learnErrors(vector_of_dereps)` aggregates the transitions.
Inference remains independent for each sample. With FASTQ paths, `learnErrors`
uses by default at most 100 million bases, split across the files, and
parallelises the raw recount across samples. Use `nbases=typemax(Int)` to
consume all the files.

On nodes with many cores, samples run as dynamic tasks: short ones do not leave
threads idle while the heavier samples finish. The alignment cache speeds up
the `learnErrors` iterations; when the RAM assigned to the job matters more
than time, use `cacheAlignments=false`. DADA2.jl uses the threads of a single
Julia process and does not communicate across nodes directly. On multiple
nodes, learn the shared model once, then distribute the per-sample inference
with that model through the job scheduler.

Qualities must be covered by the model. Dereplicated sequences must contain
only A/C/G/T; `maxN>0` in the filter can keep N, but inference requires them to
be removed. FASTA references support IUPAC codes: k-mers spanning an ambiguous
base are excluded from training. Species matches are literal, with no expansion
of ambiguity codes. Multi-line FASTQ and bzip2 are not supported.

The sequence APIs accept `LongDNA{2}`, `LongDNA{4}` and the BioSequences
`dna"..."` literals, as well as strings. This applies to alignment, FASTQ
records, chimera tests, taxonomic classification, species assignment and FASTA
writing. `readFasta(path; bio=true)` returns compact `LongDNA{4}` sequences;
the default remains `String` for compatibility. Inference and dereplication
convert the inputs once and keep the internal kernels on bytes and packed
k-mers.

To process many samples, train the indices once:

```julia
taxonomy_model = trainTaxonomy("training.fa.gz"; k=8)
taxonomy = assignTaxonomy(table.sequences, taxonomy_model; seed=42)
species_model = trainSpecies("species.fa.gz")
species = addSpecies(table.sequences, species_model)
```

## Optimisations and scalability limits

- FASTQ/gzip processed one read at a time; the whole file is never collected.
- Dereplication with a dictionary and in-place accumulation; `keepmap=false`
  drops the O(number of reads) vector when paired-end is not needed.
- K-mers encoded as UInt32 and stored as sorted O(length) multisets, without
  dense vectors of size 4^k for each sequence.
- Needleman–Wunsch with two score rows and in-band traceback:
  O(L × band) time/space for similar sequences. An upper bound on the excluded
  paths widens the band when needed; the worst case remains O(L²).
- Comparisons stored once per discovered ASV centre; O(U × C) memory, with U
  unique sequences and C ASVs, instead of an initial matrix of all pairs.
  In samples where C grows up to U, the cost remains quadratic.
- Parallelism over comparisons or over samples, avoiding nested threads.
- Log-scale probabilities, Poisson tails via Distributions and `expm1` to avoid
  numerical cancellation for small probabilities.
- Taxonomy with an inverted index of the k-mers present, instead of a dense
  4^k × taxa matrix. The database is still loaded into memory during training.
- Merging computed per distinct ASV pair rather than for every single read.

Run `julia --project=. --threads=4 benchmark/benchmarks.jl` to measure times
and allocations on your own machine. Local results and the limits of the
comparison are described in [docs/validation.md](docs/validation.md).

## Differences from DADA2 R

The Julia alignment is global with a linear penalty; upstream also uses
ends-free alignment, gapless shortcuts and homopolymer penalties.
The Julia local regression is a direct quadratic LOESS; R normally uses the
interpolated surface of `stats::loess`, so results may differ.
The comparisons do not implement the upstream greedy locks and sparse
comparison selection. Priors are handled with individual thresholds, without
reproducing all the budding priorities of the C++ code.

Chimera removal does not handle indels or reconstructions with mismatches; the
consensus uses a direct fraction of the samples containing the ASV. Taxonomy
implements its own classifier and has not been validated against the upstream
RDP classifier; `addSpecies` checks exact matches of whole sequences or
fragments.

Still to be ported: pooling/pseudo-pooling of inference, primer handling,
PhiX screening, PacBio mode, plotting functions and the other R helpers.
These differences may change ASVs and abundances on real data, even though the
verified fixture matches when the error model is shared.

## Reproducing the upstream comparison

```powershell
julia --project=. --threads=4 test/upstream_smoke.jl
Rscript validation/reference.R
julia --project=. --threads=4 validation/compare.jl
```

The comparison requires DADA2 installed in R. The R script generates the
reference error model and counts; Julia uses the same model to isolate the
comparison of inference, merging and chimeras. It does not demonstrate
equivalence of error learning or validity on other datasets.

The two required FASTQ files are included in `test/fixtures`, with provenance
and hashes. The `Pkg.test()` suite also includes this check of the complete
pipeline.

## Publishing the project on GitHub

The archive `dist/DADA2.jl-0.3.0.zip` contains the distributable project.
Extract it and use the `DADA2.jl` directory as the repository root. It includes
`.github`, tests, fixtures, examples, licence and scripts; it excludes installed
dependencies, the upstream checkout and locally generated results. The
`.sha256` file lets you verify the integrity of the archive.

Create an empty repository on GitHub, then from the extracted directory:

```sh
git init -b main
git add .
git commit -m "Release DADA2.jl 0.3.0"
git remote add origin https://github.com/alexology/DADA2.jl.git
git push -u origin main
```

Replace `YOUR_ACCOUNT` with the chosen account. After the push, GitHub Actions
runs the tests on Windows/Linux and Julia 1.10/stable. The manual
**Source archive** workflow builds a downloadable archive; it does not publish
releases. To regenerate it locally on Windows:
`powershell -NoProfile -ExecutionPolicy Bypass -File scripts/package.ps1`.
The option applies to that process only and does not change the system policy.
See [CONTRIBUTING.md](CONTRIBUTING.md) and [CHANGELOG.md](CHANGELOG.md).

Once published, other users can install the project with
`Pkg.add(url="https://github.com/alexology/DADA2.jl.git")`.
The package is not yet registered in the Julia General registry.

Licensed under LGPL-3.0-only; references and provenance in [NOTICE.md](NOTICE.md).
