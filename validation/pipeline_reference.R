# Reference outputs of R dada2 for test/fixtures/pipeline: pooled and
# pseudo-pooled inference, mergePairs and chimera detection.
# Run from the repository root: Rscript validation/pipeline_reference.R
suppressPackageStartupMessages(library(dada2))
set.seed(20261104)
out <- file.path("test", "fixtures", "pipeline")
unlink(out, recursive=TRUE)
dir.create(out, recursive=TRUE)
nt <- c("A","C","G","T")
rc <- function(s) dada2:::rc(s)
mutate <- function(x, nsub, nindel=0) {
  for(k in seq_len(nsub)) { i <- sample(length(x), 1); x[i] <- sample(setdiff(nt, x[i]), 1) }
  for(k in seq_len(nindel)) {
    i <- sample(2:(length(x)-1), 1)
    if(runif(1) < 0.5) x <- x[-i] else x <- append(x, sample(nt, 1), i)
  }
  x
}
phred <- function(maxq) {
  e <- sapply(0:maxq, function(q) min(0.75, 10^(-q/10)))
  m <- matrix(0, 16, maxq+1)
  for(f in 0:3) for(t in 0:3) m[4*f+t+1,] <- if(f == t) 1-e else e/3
  m
}
err <- phred(40)
simulate <- function(truths, abund, path, from=1, to=NULL) {
  reads <- character(0); quals <- character(0)
  for(t in seq_along(truths)) for(r in seq_len(abund[t])) {
    x <- truths[[t]]
    x <- x[from:(if(is.null(to)) length(x) else min(to, length(x)))]
    q <- pmin(40, pmax(2, round(rnorm(length(x), 30, 7) - seq_along(x)/length(x)*6)))
    e <- runif(length(x)) < 10^(-q/10)
    x[e] <- sapply(x[e], function(b) sample(setdiff(nt, b), 1))
    reads <- c(reads, paste(x, collapse="")); quals <- c(quals, rawToChar(as.raw(q+33)))
  }
  writeLines(paste0("@r", seq_along(reads), "\n", reads, "\n+\n", quals), path)
}
# Dereplicate, keep integer qualities (all the C++ core uses) and write the uniques.
derep_written <- function(fastq, path) {
  d <- derepFastq(fastq)
  d$quals <- floor(d$quals + 0.5)
  q <- apply(d$quals, 1, function(v) rawToChar(as.raw(v[!is.na(v)] + 33)))
  writeLines(paste(names(d$uniques), d$uniques, q, sep="\t"), path)
  d
}
write_dada <- function(res, prefix) {
  write.table(data.frame(sequence=names(res$denoised), abundance=unname(res$denoised)),
              paste0(prefix, "_denoised.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  writeLines(as.character(res$map), paste0(prefix, "_map.txt"))
  cl <- res$clustering
  cl$sequence <- as.character(cl$sequence)
  for(col in c("pval", "birth_pval", "birth_fold", "birth_qave"))
    cl[[col]] <- ifelse(is.na(cl[[col]]), "NA", sprintf("%.17g", cl[[col]]))
  write.table(cl, paste0(prefix, "_clustering.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}

# --- Pooling: four samples sharing rare variants -----------------------------
tmp <- tempfile(); dir.create(tmp)
base <- sample(nt, 120, TRUE)
truths <- list(base)
for(t in 2:10) truths[[t]] <- mutate(truths[[sample(t-1, 1)]], sample(c(1, 1, 2, 3), 1))
dereps <- list()
for(s in 1:4) {
  abund <- round(10^runif(10, 0, 2.3))
  abund[sample(10, 3)] <- sample(1:3, 3, TRUE)  # rare in this sample
  fq <- file.path(tmp, sprintf("pool%d.fastq", s))
  simulate(truths, abund, fq)
  dereps[[s]] <- derep_written(fq, file.path(out, sprintf("pool%d_uniques.tsv", s)))
}
write.table(matrix(sprintf("%.17g", err), 16), file.path(out, "err.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
for(mode in c("pooled", "pseudo")) {
  res <- dada(dereps, err=err, pool=if(mode == "pooled") TRUE else "pseudo",
              multithread=FALSE, verbose=FALSE)
  for(s in 1:4) write_dada(res[[s]], file.path(out, sprintf("%s%d", mode, s)))
}

# --- mergePairs ---------------------------------------------------------------
amplicon <- sample(nt, 230, TRUE)
truths <- list(amplicon)
for(t in 2:6) truths[[t]] <- mutate(truths[[sample(t-1, 1)]], sample(1:3, 1), if(t == 6) 1 else 0)
short <- c(sample(nt, 40, TRUE), amplicon[60:140])  # shorter than the reads: overhangs
truths[[7]] <- short
abund <- c(400, 150, 90, 60, 40, 30, 50)
fqF <- file.path(tmp, "F.fastq"); fqR <- file.path(tmp, "R.fastq")
set.seed(7); simulate(truths, abund, fqF, 1, 150)
# Reverse reads: reverse complement of the last 150 bases (whole short amplicon).
rtruths <- lapply(truths, function(x) strsplit(rc(paste(x[max(1, length(x)-149):length(x)], collapse="")), "")[[1]])
set.seed(8); simulate(rtruths, abund, fqR)
dF <- derep_written(fqF, file.path(out, "mergeF_uniques.tsv"))
dR <- derep_written(fqR, file.path(out, "mergeR_uniques.tsv"))
writeLines(as.character(dF$map), file.path(out, "mergeF_readmap.txt"))
writeLines(as.character(dR$map), file.path(out, "mergeR_readmap.txt"))
fF <- dada(dF, err=err, multithread=FALSE, verbose=FALSE)
fR <- dada(dR, err=err, multithread=FALSE, verbose=FALSE)
write_dada(fF, file.path(out, "mergeF")); write_dada(fR, file.path(out, "mergeR"))
settings <- list(default=list(), mismatch=list(maxMismatch=2, minOverlap=20),
                 trim=list(trimOverhang=TRUE), concat=list(justConcatenate=TRUE))
for(name in names(settings)) {
  m <- do.call(mergePairs, c(list(fF, dF, fR, dR, returnRejects=TRUE), settings[[name]]))
  m$sequence <- as.character(unlist(m$sequence))
  m$forward <- as.character(fF$clustering$sequence[m$forward])
  m$reverse <- as.character(fR$clustering$sequence[m$reverse])
  write.table(m[, c("sequence","abundance","forward","reverse","nmatch","nmismatch","nindel","prefer","accept")],
              file.path(out, paste0("merge_", name, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
}

# --- Chimeras -----------------------------------------------------------------
parents <- lapply(1:8, function(i) sample(nt, 150, TRUE))
for(i in 2:8) parents[[i]] <- mutate(parents[[1]], sample(8:20, 1))  # related parents
seqs <- sapply(parents, paste, collapse="")
chims <- character(0)
for(k in 1:30) {
  ij <- sample(8, 2); cut <- sample(30:120, 1)
  x <- c(parents[[ij[1]]][1:cut], parents[[ij[2]]][(cut+1):150])
  kind <- k %% 5
  if(kind == 1) x <- mutate(x, 1)                      # one-off chimera
  if(kind == 2) x <- mutate(x, 0, 1)                   # indel
  if(kind == 3) x <- c(sample(nt, 3, TRUE), x[1:147])  # shifted
  if(kind == 4) x <- mutate(x, 3)                      # not a chimera
  chims <- c(chims, paste(x, collapse=""))
}
allseqs <- unique(c(seqs, chims))
nsam <- 6
tab <- matrix(0L, nsam, length(allseqs), dimnames=list(paste0("s", 1:nsam), allseqs))
for(s in 1:nsam) {
  tab[s, 1:8] <- as.integer(round(10^runif(8, 1.5, 3)))
  present <- sample(9:length(allseqs), 18)
  tab[s, present] <- as.integer(round(10^runif(18, 0, 1.8)))
}
tab <- tab[, colSums(tab) > 0]
write.table(data.frame(sequence=colnames(tab), t(tab), check.names=FALSE),
            file.path(out, "chimera_table.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
flags <- data.frame(
  consensus = isBimeraDenovoTable(tab, multithread=FALSE),
  consensus_oneoff = isBimeraDenovoTable(tab, allowOneOff=TRUE, multithread=FALSE),
  consensus_loose = isBimeraDenovoTable(tab, minSampleFraction=0.5, ignoreNNegatives=2,
                                        minFoldParentOverAbundance=1, minParentAbundance=4,
                                        maxShift=4, multithread=FALSE),
  pooled = isBimeraDenovo(tab, multithread=FALSE),
  pooled_oneoff = isBimeraDenovo(tab, allowOneOff=TRUE, minOneOffParentDistance=2, multithread=FALSE))
write.table(flags, file.path(out, "chimera_flags.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
persample <- removeBimeraDenovo(tab, method="per-sample")
write.table(data.frame(sequence=colnames(persample), t(persample), check.names=FALSE),
            file.path(out, "chimera_persample.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
writeLines(c(paste("R", getRversion()), paste("dada2", packageVersion("dada2"))),
           file.path(out, "versions.txt"))
