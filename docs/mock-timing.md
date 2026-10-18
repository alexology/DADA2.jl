# Timing on the mock community: Julia versus DADA2 R

**The complete Julia pipeline is slower than R on this sample.**
The local comparison of 10 October 2026 uses Julia 1.12.2, DADA2.jl 0.3.0,
R 4.5.0 and DADA2 R 1.36.0, on Windows with an Intel Core i7-9750H, 6 cores and
12 logical processors. The processes are run sequentially.

Public sample Mock_S280, 4,779 read pairs. Filtering at 240/160 bases,
`maxEE=2`, `maxN=0`, `truncQ=2`, `minLen=20`, PhiX disabled in R; gzip in
both. 4,314 pairs are retained. Merging requires at least 12 bases and zero
mismatches, selects amplicons of 250–256 bases and uses pooled chimera removal
with `minFoldParentOverAbundance=1`.

## Warm results

Medians of five repetitions, in seconds. Each process first runs one call
excluded from the medians, then the five measured repetitions.

| Mode | Threads requested | R | Julia |
| --- | ---: | ---: | ---: |
| Complete pipeline, mean qualities | 1 | 2.140 | 2.571 |
| Complete pipeline, mean qualities | 4 | 1.720 | 2.359 |
| Same pre-estimated error model | 1 | 1.160 | 0.927 |
| Same pre-estimated error model | 4 | 1.190 | 0.644 |
| Complete pipeline, original qualities | 1 | — | 3.179 |
| Complete pipeline, original qualities | 4 | — | 3.083 |

With four threads Julia takes **about 37% more time** for the complete
pipeline on mean qualities. Range over the five runs:
R 1.710–1.800 seconds, Julia 2.325–2.402 seconds.
The Julia variant that learns from the original qualities is reported
separately: R does not implement this mode in the method compared here.

With the errors already provided, Julia is ahead on the total: about **1.85×**
with four threads. This total covers filtering, dereplication, inference,
merging and chimeras. It does not mean that Julia's ASV inference is faster,
nor that all the functionality of the two packages is equivalent.

In **all 60 runs**, the final sequences and counts are identical to the R
reference result: **20 ASVs and 4,269 final reads**. The identity of every
sequence and its count was verified, not just the totals.

## Where the time goes

Complete pipeline on mean qualities, four threads; medians per stage:

| Stage | R, seconds | Julia, seconds |
| --- | ---: | ---: |
| Filtering | 0.450 | 0.244 |
| Dereplication | 0.100 | 0.110 |
| Forward + reverse error learning | 0.580 | 1.722 |
| Forward + reverse inference | 0.170 | 0.289 |
| Pair merging | 0.030 | 0.000218 |
| Table, length selection and chimeras | 0.390 | 0.000380 |

The medians of the individual stages do not necessarily add up to the median of
the total. The main cost in Julia is **error learning**.
In this implementation `learnErrors` parallelises inference over the vector of
samples; with the single mock sample, this path disables the internal
parallelism of the comparisons. This explains a structural limit of the
four-thread test; it does not by itself account for the whole gap with R.

The very low time for Julia chimeras must be read together with the method:
Julia checks exact recombinations between sequences of the same length, whereas
the R method has more general functionality. On this mock the result matches;
the speed advantage of this stage does not prove general equivalence.

## First run

Time from the first statement of the script to the first result of the
complete pipeline, including package loading and Julia compilation:

| Threads | R | Julia |
| --- | ---: | ---: |
| 1 | 9.52 s | 14.07 s |
| 4 | 9.36 s | 13.91 s |

The timer does not include starting the executable, installing dependencies or
a first precompilation of the package. The packages were already installed and
precompiled. The initial calls to the other modes are recorded as
`first_call`, but they share the same, already running process: they are not
independent start-ups.

## Settings and limitations

In the complete comparison on mean qualities, both packages reuse the
dereplicated objects during learning, with at most 10 iterations and
`OMEGA_C=0` during learning; final inference uses the default
`OMEGA_C=1e-40`. Julia performs three updates per direction; the counter of the
R estimation function records four calls, including the initialisation stage.
Initialisation, convergence criterion and smoothing are not identical: Julia
uses direct LOESS and R an interpolated surface.
The comparison with shared errors loads the same two previously saved R models.

R uses `multithread=FALSE` or `multithread=4` for learning/inference;
filtering of the single sample uses `multithread=FALSE` and enables OMP only in
the four-thread configuration. Julia uses processes with `--threads=1` and
`--threads=4`; BLAS is set to one thread. The requested number does not
guarantee that every stage uses all the threads.

The timers measure wall time, including reading and writing the filtered FASTQ
files. Creating/removing temporary directories, forced garbage collection and
writing the benchmark results are kept outside the pipeline time.
The final benchmark results are exported after the measurement. Peak memory
and performance on millions of reads are not measured. The operating system
cache is generally warm; the order of the processes is fixed and these are
local measurements with only five repetitions, without statistical tests.

R emits deprecation warnings for `S4Vectors::anyMissing`, kept in the files
`benchmark/mock_timing_R_*_warnings.txt`; they are not execution errors.

## Data and reproduction

- [All measurements](../benchmark/mock_timing_results.tsv)
- [Summary per stage](../benchmark/mock_timing_summary.tsv)
- [Environment and input/code hashes](../benchmark/mock_timing_metadata.json)

With the public inputs and the R models from the validation already present:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/mock-timing.ps1 -RscriptPath "C:\Program Files\R\R-4.5.0\bin\Rscript.exe"
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" --vanilla benchmark/mock_timing_summary.R
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/mock-timing-metadata.ps1
```

The counts of each repetition are in `validation/generated/mock_timing/`,
excluded from the source archives. The public data and the scientific settings
are described in [scientific-protocol.md](scientific-protocol.md).
