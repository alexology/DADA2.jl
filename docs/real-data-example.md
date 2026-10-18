# Trying DADA2.jl on real data

The script [examples/real_data.jl](../examples/real_data.jl) reads single-end
or paired-end FASTQ files, plain or gzip, with Phred+33 qualities. Run the
commands from the package directory. Quote paths that contain spaces.

## Installation and first run

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. --threads=4 examples/real_data.jl --help
julia --project=. --threads=4 examples/real_data.jl --samples examples/samples.tsv --out results_test --trunc-f 240 --trunc-r 160
```

The first run uses the included FASTQ files. The output directory must be new
or empty. To repeat a run, choose another directory name.

## A single real sample

Paired-end:

```powershell
julia --project=. --threads=4 examples/real_data.jl --forward "D:/data/sample_R1.fastq.gz" --reverse "D:/data/sample_R2.fastq.gz" --sample sample --out results_sample --trunc-f 240 --trunc-r 160 --max-ee-f 2 --max-ee-r 2
```

The lengths 240 and 160 are an example used in the included tests, not values
suitable for every dataset. Choose them according to the read length and the
quality profile. For an amplicon of length L after primer removal, the sum of
the retained lengths of the two mates must allow at least `--min-overlap`
bases of overlap. Take into account `trimLeft`, the length variability of the
amplicon and the quality cuts.
The filter discards reads that do not reach `truncLen`.

Single-end:

```powershell
julia --project=. --threads=4 examples/real_data.jl --forward "D:/data/sample.fastq.gz" --sample sample --out results_single --trunc-f 200 --max-ee-f 2
```

When `--trunc-f`/`--trunc-r` are omitted, the filter keeps variable lengths,
subject to quality truncation and the minimum length. The script does not
automatically remove primers, adapters or PhiX; use already prepared input.
`--trim-left-f` and `--trim-left-r` remove a fixed number of bases; they do not
search for primer sequences.

## Multiple samples

Prepare a TSV file, with real tab characters between the three columns:

```text
sample	forward	reverse
sample_A	data/A_R1.fastq.gz	data/A_R2.fastq.gz
sample_B	data/B_R1.fastq.gz	data/B_R2.fastq.gz
```

Relative paths are resolved against the TSV directory. Sample names must be
unique and contain ASCII letters, digits, `_`, `-` or a dot, starting with a
letter or a digit. For single-end data, leave the third column empty, keeping
the tab after the forward path. Use either all single-end or all paired-end
samples in the same run.

```powershell
julia --project=. --threads=4 examples/real_data.jl --samples samples.tsv --out results_batch --trunc-f 240 --trunc-r 160
```

The script learns one shared forward model and, when needed, a separate reverse
model, using the original FASTQ qualities. ASV inference stays separate for
each sample. Group samples from a homogeneous sequencing process. Samples with
zero filtered reads remain in the final table with zero counts and do not take
part in learning. If all samples are empty, the script stops and leaves
`filtering.tsv` for diagnosis.

## Parameters and results

`--help` lists all options. A range of ASV lengths can be selected with
`--min-asv-length` and `--max-asv-length`. For paired-end data this range
applies to the merged sequences, for single-end data to the sequences retained
after filtering. `--min-overlap` defaults to 12 and merging requires zero
mismatches. `--max-iter` defaults to 10. Chimera removal uses consensus by
default; `--chimera-method pooled` is available.

The output directory contains:

| File | Content |
| --- | --- |
| `asvs.fasta` | Final sequences, after the length and chimera checks |
| `asv_counts.tsv` | Rows = samples, columns = ASVs, values = counts |
| `asvs_all.fasta`, `asv_counts_all.tsv` | Results before the length and chimera checks |
| `asv_summary.tsv` | Length, total count and kept/length/chimera status |
| `read_tracking.tsv` | Counts after filtering, inference, merging, length and chimeras |
| `filtering.tsv` | Filtering summary, written even if the analysis stops after that step |
| `errors_forward.tsv`, `errors_reverse.tsv` | Learned models; reverse only for paired-end |
| `error_support_*.tsv`, `learning_*.tsv` | Exposures per base/quality and convergence history |
| `run.toml` | Versions, parameters, inputs and convergence status |
| `filtered/` | Filtered FASTQ files, useful to repeat or compare the analysis |

ASV identifiers stay consistent across all the files of the same run; after
exclusions some numbers may be missing. They are assigned for that run and are
not global identifiers across different analyses.

For paired-end data, `reads_in`, `reads_filtered`, `merged` and `nonchimeric`
count pairs/merged sequences, not the sum of forward and reverse reads.
A sharp drop in `merged` may be due to insufficient overlap.
Non-convergence is reported and recorded: the files produced should be treated
as preliminary results. `supported=true` means at least one observation, not
that the estimate is precise.

The pipeline limitations are those described in the README: merging without
indels and exact chimera checks between sequences of the same length.
The script produces ASVs and abundances; it does not assign taxonomy
automatically.

Verified locally with Julia 1.12.2: paired-end example run directly and from a
TSV (7 ASVs, 938 final reads), single-end (10 ASVs, 1,083 final reads) and a
batch with two active samples and one empty sample. Check the consistency of
the FASTA/TSV identifiers, the count sums and the retention of reads. The batch
learns a different model from the single sample, so it may produce different
counts even for the same input.
