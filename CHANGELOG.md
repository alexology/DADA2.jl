# Changelog

## Unreleased

- BioSequences 3.5 interoperability across DNA alignment, FASTQ records,
  chimera detection, taxonomy/species queries and FASTA writing.
- `readFasta(path; bio=true)` returns compact `LongDNA{4}` values while the
  default string API remains compatible.
- FASTQ error learning now has a 100-million-base default training budget,
  split deterministically across sample prefixes; `nbases=typemax(Int)` keeps
  the previous all-data behavior.
- Raw-quality transition recounting runs concurrently across samples with
  thread-local matrices and a deterministic reduction. Results and diagnostics
  expose the selected read and base counts.
- Error learning defaults to `omegaC=0`, so assigned reads contribute to error
  estimation as in the upstream learning procedure; callers can override it.

## 0.3.0 — 2026-10-10

- Real-data CLI example for single/paired-end and sample TSV manifests,
  with shared error learning, read tracking, diagnostics and FASTA/TSV exports.

- Native direct quadratic LOESS with exposure weights, constant boundary
  rates and stable lower-degree fits for sparse quality support. The previous
  linear estimator remains selectable with `smoother=:linear`.
- FASTQ-based `learnErrors` counts transitions at original quality scores
  by default (`qualityMode=:raw`), streaming each iteration. Derep inputs
  retain mean-quality counting; FASTQ callers may select `qualityMode=:mean`.
- Error-learning results expose transition counts, inference convergence and
  quality mode; `errorDiagnostics` reports training exposures and substitutions.
- Independent R smoothing fixtures and scientific calibration comparisons,
  including shared models across five simulated samples and experimental mock.

## 0.2.1 — 2026-10-10

- Reject reads shortened below `truncLen` by quality trimming, before `trimLeft`,
  matching DADA2 R. Paired filtering rejects the entire pair if either read fails.
- Add regression tests for quality cutoff boundaries and left trimming.

## 0.2.0 — 2026-10-10

- Added preliminary scientific validation: public experimental mock community,
  30 simulated conditions, R comparison, raw metrics, protocol and figures.
  No denoiser changes were made during validation; known failures are documented.

- Reusable sparse taxonomy and species indexes via `trainTaxonomy` and `trainSpecies`.
- IUPAC reference FASTA and reverse complements; k-mers skip ambiguous positions.
- Exact species matching against fragments of longer references, using a seed index.
- Linear selection of distinct chimera parents instead of a quadratic parent-pair search.
- Constant-time center membership, cached log abundances, and no alignment of centers to themselves.
- Robust error-learning convergence calculation when transition probabilities are zero.
- Bundled upstream paired-end fixtures, integration tests and offline source archives.
- GitHub CI, source archive workflow, contribution instructions and release packaging.

## 0.1.0 — 2026-10-10

- Initial Julia pipeline: FASTQ/gzip, filtering, dereplication, ASV inference,
  error learning, merging, chimera detection, taxonomy and FASTA IO.
- Shared-error-model parity on the upstream sam1 fixture with DADA2 R 1.36.0.
