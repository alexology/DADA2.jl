# Reference outputs of the DADA2 C++ core for test/fixtures/core.
# Run from the repository root: Rscript validation/core_reference.R
# Simulated uniques are denoised by R dada2 with fixed error matrices and one
# option set per case; random pairs are aligned with the upstream aligners.
suppressPackageStartupMessages(library(dada2))
set.seed(20261102)
out <- file.path("test", "fixtures", "core")
unlink(out, recursive=TRUE)
dir.create(out, recursive=TRUE)
nt <- c("A","C","G","T")

mutate <- function(x, nsub, nindel) {
  for(k in seq_len(nsub)) { i <- sample(length(x), 1); x[i] <- sample(setdiff(nt, x[i]), 1) }
  for(k in seq_len(nindel)) {
    i <- sample(length(x), 1)
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

# Alignments: vectorized/classic, ends-free/global, homopolymer gaps, bands.
rows <- list()
for(k in 1:400) {
  len <- sample(c(6:40, 100:180), 1)
  s1 <- sample(nt, len, TRUE)
  kind <- sample(3, 1)
  s2 <- if(kind == 1) sample(nt, sample(6:180, 1), TRUE) else
        if(kind == 2) mutate(s1, sample(0:6, 1), sample(0:4, 1)) else
        c(rep(sample(nt, 1), sample(3:6, 1)), mutate(s1, 1, 1))
  s1 <- paste(s1, collapse=""); s2 <- paste(s2, collapse="")
  if(runif(1) < 0.5) { tmp <- s1; s1 <- s2; s2 <- tmp }
  band <- sample(c(-1, 1, 2, 3, 4, 5, 8, 16, 16, 32), 1)
  endsfree <- runif(1) < 0.7
  vec <- runif(1) < 0.6
  homo <- if(!vec && endsfree && runif(1) < 0.5) -1L else -8L
  al <- nwalign(s1, s2, match=5L, mismatch=-4L, gap=-8L, homo_gap=homo, band=band,
                endsfree=endsfree, vec=vec)
  rows[[k]] <- data.frame(s1, s2, band, endsfree, vec, homo, al0=al[1], al1=al[2])
}
write.table(do.call(rbind, rows), file.path(out, "nwalign.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Partitions: one option set per case.
options <- list(list(), list(OMEGA_C=0), list(GREEDY=FALSE), list(GAPLESS=FALSE),
                list(VECTORIZED_ALIGNMENT=FALSE), list(HOMOPOLYMER_GAP_PENALTY=-1),
                list(BAND_SIZE=-1), list(BAND_SIZE=0), list(USE_KMERS=FALSE),
                list(DETECT_SINGLETONS=TRUE), list(MIN_FOLD=2, MIN_HAMMING=2, MIN_ABUNDANCE=3),
                list(OMEGA_A=1e-10, KDIST_CUTOFF=0.2), list(MAX_CLUST=3), list(PRIORS=TRUE))
maxq <- 40
for(k in seq_along(options)) {
  dir <- file.path(out, sprintf("case%02d", k)); dir.create(dir)
  len <- sample(c(60, 90, 120), 1)
  truths <- list(sample(nt, len, TRUE))
  ntrue <- sample(4:9, 1)
  for(t in 2:ntrue) truths[[t]] <- mutate(truths[[sample(t-1, 1)]], sample(c(1, 1, 2, 4), 1),
                                          if(runif(1) < 0.3) 1 else 0)
  abund <- round(10^runif(ntrue, 0.3, 2.7))
  reads <- character(0); quals <- character(0)
  for(t in seq_len(ntrue)) for(r in seq_len(abund[t])) {
    x <- truths[[t]]
    q <- pmin(maxq, pmax(2, round(rnorm(length(x), 30, 7) - seq_along(x)/length(x)*8)))
    err <- runif(length(x)) < 10^(-q/10)
    x[err] <- sapply(x[err], function(b) sample(setdiff(nt, b), 1))
    if(runif(1) < 0.01) x <- mutate(x, 0, 1)
    q <- q[seq_along(x)]; q[is.na(q)] <- 30
    reads <- c(reads, paste(x, collapse="")); quals <- c(quals, rawToChar(as.raw(q+33)))
  }
  fq <- tempfile(fileext=".fastq")
  writeLines(paste0("@r", seq_along(reads), "\n", reads, "\n+\n", quals), fq)
  drp <- derepFastq(fq)
  # The C++ core only uses qualities rounded half away from zero; store those.
  drp$quals <- floor(drp$quals + 0.5)
  err <- phred(maxq)
  if(k %% 2 == 0) {
    err <- err * matrix(runif(16*(maxq+1), 0.3, 3), 16)
    for(f in 0:3) {
      self <- 4*f+f+1; others <- setdiff(4*f+(1:4), self)
      err[others,] <- pmin(err[others,], 0.2); err[self,] <- 1 - colSums(err[others,,drop=FALSE])
    }
  }
  opt <- options[[k]]
  priors <- character(0)
  if(!is.null(opt$PRIORS)) { opt$PRIORS <- NULL; priors <- sapply(truths, paste, collapse="") }
  res <- do.call(dada, c(list(drp, err=err, multithread=FALSE, verbose=FALSE, priors=priors), opt))
  q <- apply(drp$quals, 1, function(v) rawToChar(as.raw(v[!is.na(v)] + 33)))
  writeLines(paste(names(drp$uniques), drp$uniques, q, sep="\t"), file.path(dir, "uniques.tsv"))
  write.table(matrix(sprintf("%.17g", err), 16), file.path(dir, "err.tsv"), sep="\t",
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  writeLines(as.character(sapply(seq_along(opt), function(i)
    paste(names(opt)[i], format(opt[[i]], digits=17), sep="="))), file.path(dir, "options.txt"))
  if(length(priors)) writeLines(priors, file.path(dir, "priors.txt"))
  write.table(data.frame(sequence=names(res$denoised), abundance=unname(res$denoised)),
              file.path(dir, "denoised.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
  writeLines(as.character(res$map), file.path(dir, "map.txt"))
  write.table(res$trans, file.path(dir, "trans.tsv"), sep="\t", row.names=FALSE, col.names=FALSE)
  cl <- res$clustering
  cl$sequence <- as.character(cl$sequence)
  for(col in c("pval", "birth_pval", "birth_fold", "birth_qave"))
    cl[[col]] <- ifelse(is.na(cl[[col]]), "NA", sprintf("%.17g", cl[[col]]))
  write.table(cl, file.path(dir, "clustering.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}
writeLines(c(paste("R", getRversion()), paste("dada2", packageVersion("dada2"))),
           file.path(out, "versions.txt"))
