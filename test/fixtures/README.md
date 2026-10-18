# Upstream fixtures

`sam1F.fastq.gz` and `sam1R.fastq.gz` are copied without modification from
https://github.com/benjjneb/dada2/tree/b61fe2e90115be08e78a01a16af89fd21402ca42/inst/extdata.
The original DADA2 license is preserved in the root LICENSE; see NOTICE.md
for attribution and the upstream licensing metadata discrepancy.

SHA-256:

```
sam1F.fastq.gz f670b674d9a066314760a7ff0b3850264602b0050d19692438c3fe179efa6544
sam1R.fastq.gz ee11adfb81fb0cb4f972b92bedb91a9e94ac694d2900ade4111421ad220d7372
```

These fixtures contain 1,500 paired reads and allow integration tests to run
offline after Julia dependencies have been installed.

## Error smoothing fixtures

`error_transitions.tsv` is synthetic data with curved log error rates,
unequal exposures and unobserved quality tails, generated independently
by `validation/error_smoothing_fixture.jl`.
`error_loess_direct.tsv` and `error_loess_r.tsv` are reference predictions
from R's direct LOESS and DADA2's default interpolated LOESS, respectively.
Regenerate with `validation/error_smoothing_reference.R`; software versions
are in `error_loess_versions.txt`. Tests use the direct reference and need
no R installation at runtime.
