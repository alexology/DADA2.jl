# Provenance

Julia implementation based on the DADA2 model and pipeline by Benjamin
Callahan and collaborators: https://github.com/benjjneb/dada2.
Sources consulted at commit `b61fe2e90115be08e78a01a16af89fd21402ca42`
(DESCRIPTION: version 1.36.0). The Julia code was written in this project;
the inference core, aligners, error-model smoothing, pair merging, chimera
and phiX screening are translations of the upstream R/C++ code (and of R's
`loess` Fortran for the one-predictor case). R's Mathlib code is not used.
`data/phix_genome.fa` is copied unmodified from the upstream
`inst/extdata/phix_genome.fa` (SHA-256
19beaeb235899980ec5019a475e25fc41dd6ea98a7bf7afb51dace28b73cc19a).
The fixtures `test/fixtures/sam1F.fastq.gz` and `sam1R.fastq.gz` are reproduced
unmodified from that upstream commit, with provenance and checksums in the
fixtures README.

The LICENSE file reproduces the LGPL version 3 present in the upstream
repository. The License field of its DESCRIPTION states LGPL-2: this
discrepancy is recorded here to preserve provenance. The new code of this
project is distributed under LGPL-3.0-only.

Scientific reference: Callahan et al., “DADA2: High-resolution sample inference
from Illumina amplicon data”, Nature Methods 13, 581–583 (2016),
https://doi.org/10.1038/nmeth.3869.
