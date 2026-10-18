# DADA2.jl 0.2.1

A read could pass the initial `truncLen` check, become shorter after quality
truncation and still be kept. The filter now discards it. For paired-end data,
failure of one mate is enough to discard the whole pair.

The check compares the length before subtracting `trimLeft`:
`truncLen=6, trimLeft=2` produces four-base reads. With `truncLen=0`, variable
length remains allowed, subject to `minLen` and `maxLen`.

Local verification of 10 October 2026, Julia 1.12.2 on Windows:

- All **622 assertions** of the suite pass, including 16 new regression tests
  on early truncation, quality threshold, truncation position, initially short
  reads, `trimLeft` and failure of each mate.
- **6 additional checks** on the public sample Mock_S280: 4,779 initial pairs,
  **4,314 retained**, as in DADA2 R 1.36.0. Identifiers, order, sequences and
  qualities match for both mates; lengths are exactly 240 and 160 bases. The
  nine extra pairs kept by 0.2.0 are now removed.

Result: [mock_filter_regression.tsv](../validation/mock_filter_regression.tsv).
To repeat the comparison, with the public data and the R outputs prepared
by `scripts/scientific.ps1`:

```powershell
julia --project=. --threads=4 -e 'using Pkg; Pkg.test()'
julia --project=. validation/mock_filter_regression.jl
```

The script uses temporary files. The scientific report, the ASV counts, the
simulations and the audit of 0.2.0 remain historical results of that version.
This verification concerns filtering; it is not a new validation of the whole
pipeline, nor of every combination of filter options.
