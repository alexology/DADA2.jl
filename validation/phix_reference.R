# Reference outputs of R dada2 isPhiX for test/fixtures/phix.tsv.
# Run from the repository root: Rscript validation/phix_reference.R
suppressPackageStartupMessages(library(dada2))
set.seed(20261103)
nt <- c("A","C","G","T")
phix <- as.character(ShortRead::sread(ShortRead::readFasta(system.file("extdata", "phix_genome.fa", package="dada2"))))
g <- strsplit(phix, "")[[1]]
mutate <- function(x, n) { for(k in seq_len(n)) { i <- sample(length(x), 1); x[i] <- sample(nt, 1) }; x }
seqs <- character(0)
for(k in 1:600) {
  kind <- k %% 6
  len <- sample(c(10:40, 60, 120, 250), 1)
  x <- if(kind == 0) sample(nt, len, TRUE) else {
    start <- if(kind == 5) length(g) - sample(5:30, 1) else sample(length(g), 1)
    idx <- ((start + seq_len(len) - 2) %% length(g)) + 1  # circular
    mutate(g[idx], sample(c(0, 1, 3, 8, 20), 1))
  }
  if(kind == 3) x <- strsplit(dada2:::rc(paste(x, collapse="")), "")[[1]]
  if(kind == 4 && length(x) > 20) x[sample(length(x), 2)] <- "N"
  if(kind == 2 && len >= 40) x <- c(sample(nt, 20, TRUE), x[1:min(40, length(x))], sample(nt, 20, TRUE))
  seqs <- c(seqs, paste(x, collapse=""))
}
out <- data.frame(sequence=seqs,
  default=isPhiX(seqs),
  w12_m3=isPhiX(seqs, wordSize=12, minMatches=3),
  overlapping=isPhiX(seqs, nonOverlapping=FALSE, minMatches=3),
  w20_m1=isPhiX(seqs, wordSize=20, minMatches=1))
write.table(out, file.path("test", "fixtures", "phix.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
print(colSums(out[, -1]))
