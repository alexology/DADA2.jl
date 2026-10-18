# Provenance

Julia implementation based on the DADA2 model and pipeline by Benjamin
Callahan and collaborators: https://github.com/benjjneb/dada2.
Sources consulted at commit `b61fe2e90115be08e78a01a16af89fd21402ca42`
(DESCRIPTION: version 1.36.0). The Julia code was written in this project;
it does not incorporate the C++ kernels or R's Mathlib code.
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
