# Scientific validation protocol

Date: 10 October 2026. Code DADA2.jl 0.2.0; R reference DADA2 1.36.0.
Parameters were defined in the scripts before the results were analysed; no
change to the denoiser or the thresholds will be made to improve these results.

## Simulated experiments

Six scenarios × five seeds (11, 23, 37, 51, 71), 30 simulated conditions:
no errors, calibrated errors, rare SNVs, miscalibrated qualities, indels and
chimeras. The reference consists of five 150-base sequences, including two
variants one substitution away from an abundant sequence.
The sequences are generated at random and the counts are known before errors
are introduced. The five seeds define five independent communities; the same
seed is reused across scenarios, so the 30 cases are not 30 independent
replicates.

Abundances are normally 3000/1000/500/100/20; in the rare scenario
5000/1000/200/10/1. Qualities are sampled uniformly in 25:40.
Transitions have asymmetric weights 1:2:3, built independently of the
package's `phredErrors` model. In the miscalibrated case the error probability
is five times the nominal one. Random indels have probability 0.001 per
position for each of the insertion/deletion processes; a fifth ASV carries a
true insertion. In the chimera case, 200 reads of an exact recombination of the
first two parents are added.

Comparisons: Julia/R with the known model; Julia/R with independent learning;
Julia with the errors learned by R. Default inference parameters, singletons
disabled; pooled chimeras, minimum parent ratio 1. No ASV pooling across
replicates. No filtering that would change the expected biological count.

Metrics: TP/FP/FN on exact sequence identity, precision, recall, recall
excluding biological singletons, total variation distance between relative
abundances (half the sum of the absolute differences), fraction of reads
assigned to spurious sequences and Julia convergence.
Discarded reads are reported indirectly through the accepted counts; the
distance on relative abundances does not by itself measure absolute loss.

A singleton ASV contributes to full recall but not to `recall_ge2`.
This avoids holding the method to a target that the default option excludes.
Artificial chimeras are not part of the biological ground truth.
A known model does not make the simulator realistic, nor does it demonstrate
validity for correlated errors, PCR, index hopping or experimental
contamination.

## Experimental mock community

Public mothur MiSeq SOP dataset, also used in the
[DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).
Download: https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip.
Only Mock_S280 and the included HMP_MOCK.v35.fasta reference are analysed.

Common filtering: truncLen 240/160, maxEE 2/2, maxN 0, truncQ 2.
PhiX disabled in both pipelines because Julia does not implement that check.
Learning on this sample, without using the mouse samples; merging with
minOverlap 12, zero mismatches; amplicon lengths 250:256 as defined by the
tutorial; pooled chimeras, minimum parent ratio 1.

Metrics: ASVs exactly matching a fragment of the reference, non-matching ASVs,
precision with respect to the reference, fraction of non-matching reads, number
of reference records covered. The records are not assumed to all be distinct
ASVs, and recall is not computed over the 21 nominal strains: the tutorial
reports one strain as absent. The true abundances of the individual amplicons
are unknown, so count bias is not computed.

The references may be incomplete: "non-matching" is an operational definition
and does not by itself prove that the sequence is a biological false positive.
Precision with respect to the references is not equivalent to taxonomic
accuracy.

## Interpretation

These experiments constitute a preliminary quantitative validation, not a
certification of general equivalence. Five replicates per scenario make it
possible to report variability, not narrow intervals or conclusions about all
microbiomes. Negative results remain in the report and in the TSV files.
A production-grade validation requires further independent mocks, platforms,
depths, amplicon regions, negative controls and taxonomic benchmarks.
