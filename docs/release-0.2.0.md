# DADA2.jl 0.2.0

Release of the Julia package with a FASTQ → ASV pipeline, learnable error
models, paired-end merging and chimera checking. This implementation remains
experimental and retains the scientific differences described in the README.

Improvements in this release:

- taxonomic and species indices reusable across samples;
- IUPAC references and species assigned to exact fragments;
- more efficient chimera checks and handling of inference centres;
- public fixtures included, pipeline tests and GitHub workflow;
- reproducible source archive with SHA-256 checksum;
- preliminary validation on a real mock and 30 simulated conditions,
  with documented individual results and scientific limitations.

Installation after uploading to your own repository:

```julia
using Pkg
Pkg.add(url="https://github.com/YOUR_ACCOUNT/DADA2.jl.git", rev="v0.2.0")
```

To use this command, first create the `v0.2.0` tag in your own repository.
The archive and the tag do not constitute registration in Julia General.
Comparisons with R and the limits of the model are described in `docs/validation.md`.
