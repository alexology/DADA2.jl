# Local checks — 10 October 2026

0.3.0 update: **662 assertions pass**. The new error learning is verified in
[error-learning.md](error-learning.md), with an independent LOESS comparison,
original qualities, shared models and scientific metrics.
The 0.2.0 measurements reported below still refer to that version.

Environment: Windows, Julia 1.12.2, four threads; R 4.5.0 with DADA2 1.36.0.
Upstream consulted: commit `b61fe2e90115be08e78a01a16af89fd21402ca42`.
The measurements are local and not guaranteed on other machines or other data.

## Tests

Version 0.2.0: **606 assertions pass**, also with the standard command
`julia --project=. --threads=4 -e 'using Pkg; Pkg.test()'`.
Covered: FASTQ/gzip, paired-end filters and cleanup on error, dereplication
qualities and maps, 120 alignments compared with an independent dense DP,
small Poisson tails, ASV budding, OMEGA_C correction, convergence even with
zero probabilities, serial/parallel results, merging, chimeras compared with an
exhaustive check, reusable taxonomy, IUPAC references, species assignment to
fragments and FASTA.

The alignment tests check the optimal global score; they do not check that the
traceback is identical to the upstream ends-free kernels.
The CI workflow is set up for Julia 1.10 and stable, Windows and Linux;
these remote runs were not performed in the local session.

## Comparison with DADA2 R on the sam1 fixture

Filtering: `truncLen=(240,160)`, `maxEE=(2,2)`, 1,500 input pairs,
1,021 retained by both implementations.

Models estimated by DADA2 R and reused in Julia:

| Stage | R ASVs | Julia ASVs | R reads | Julia reads | Identical sequences and counts |
| --- | ---: | ---: | ---: | ---: | --- |
| Forward inference | 9 | 9 | 1,017 | 1,017 | Yes |
| Reverse inference | 8 | 8 | 1,021 | 1,021 | Yes |
| Merging and chimera removal | 7 | 7 | 965 | 965 | Yes |

The four forward reads excluded by the OMEGA_C check do not contribute to the
final abundances. Abundance-weighted Jaccard: 1.0 for both forward and reverse.

Pipeline with learning done entirely in Julia, at most ten iterations:
convergence reached for both models and both inferences, 9 forward ASVs,
8 reverse and 7 final, but **938 merged reads**, against 965 using the errors
estimated by R. The difference confirms that Julia's learning is not
equivalent to the upstream loess.

The initial comparison did not constitute a general scientific validation.
The later experiments on a real mock community and 30 simulated conditions are
described in [scientific-validation.md](scientific-validation.md), including
the loss of rare variants and the filtering differences that emerged.

## Synthetic benchmark

Measured after compilation, minimum of five repetitions with GC before each.
The alignment compares two 1,000-base sequences with one substitution.
The dense reference computes only the score, while the optimised Julia code
also produces the alignment and the position map.

| Operation | Minimum time | Allocated bytes |
| --- | ---: | ---: |
| Banded alignment | 0.304 ms | 75,516 |
| Dense DP, score only | 9.306 ms | 8,016,167 |
| Serial inference, 100 sequences of 250 bases | 42.276 ms | 9,306,383 |
| Four-thread inference, same input | 32.890 ms | 9,495,263 |
| Taxonomy training and classification | 5.086 ms | 7,719,519 |
| Classification with a pre-trained index | 0.649 ms | 469,006 |

In this case the alignment allocates about 106 times fewer bytes and is about
31 times faster than the dense Julia reference. The inference measurements show
a modest improvement with four threads: the cost of budding and of the
reassignment cycles remains serial. This is not a Julia versus DADA2 R/C++
benchmark and it does not measure the peak memory of the process.

The taxonomy comparison uses 100 references of 250 bases, classification of
10 sequences and 10 bootstrap replicates. Reusing the index is about 7.8 times
faster and allocates about 16.5 times fewer bytes than repeated training and
classification. The initial training cost does not disappear: it is paid once.

The files `benchmark/results.tsv` and `validation/generated/` can be
regenerated with the scripts in the repository. The two required FASTQ files
are included in `test/fixtures`; the local upstream checkout is excluded from
the archive.
