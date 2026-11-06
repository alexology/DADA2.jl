# Reference outputs of R learnErrors for test/fixtures/learn.
# Run from the repository root: Rscript validation/learn_reference.R
# Three small simulated samples; learnErrors on all of them and with an
# nbases budget that stops after the second file.
suppressPackageStartupMessages(library(dada2))
set.seed(20261103)
out <- file.path("test", "fixtures", "learn")
unlink(out, recursive=TRUE)
dir.create(out, recursive=TRUE)
nt <- c("A","C","G","T")
mutate <- function(x, nsub) {
  for(k in seq_len(nsub)) { i <- sample(length(x), 1); x[i] <- sample(setdiff(nt, x[i]), 1) }
  x
}
files <- character(0)
for(s in 1:3) {
  truths <- list(sample(nt, 100, TRUE))
  for(t in 2:6) truths[[t]] <- mutate(truths[[sample(t-1, 1)]], sample(1:4, 1))
  abund <- round(10^runif(6, 0.5, 2.5))
  reads <- character(0); quals <- character(0)
  for(t in seq_along(truths)) for(r in seq_len(abund[t])) {
    x <- truths[[t]]
    q <- pmin(40, pmax(2, round(rnorm(100, 32, 6) - seq_len(100)/12)))
    err <- runif(100) < 10^(-q/10)
    x[err] <- sapply(x[err], function(b) sample(setdiff(nt, b), 1))
    reads <- c(reads, paste(x, collapse="")); quals <- c(quals, rawToChar(as.raw(q+33)))
  }
  path <- file.path(out, sprintf("sample%d.fastq.gz", s))
  con <- gzfile(path, "w")
  writeLines(paste0("@s", s, "r", seq_along(reads), "\n", reads, "\n+\n", quals), con)
  close(con)
  files <- c(files, path)
}
bases <- sapply(files, function(f) { d <- derepFastq(f); sum(d$uniques * nchar(names(d$uniques))) })
budgets <- c(all=1e8, two=bases[[1]] + 1)
for(name in names(budgets)) {
  e <- learnErrors(files, nbases=budgets[[name]], multithread=FALSE, verbose=FALSE)
  write.table(matrix(sprintf("%.17g", e$err_out), 16), file.path(out, paste0(name, "_err.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(e$trans, file.path(out, paste0(name, "_trans.tsv")), sep="\t",
              row.names=FALSE, col.names=FALSE)
  writeLines(c(format(budgets[[name]], scientific=FALSE), length(e$err_in)),
             file.path(out, paste0(name, "_info.txt")))
}
# Pooled and pseudo-pooled learning on all three samples.
for(mode in c("pooled", "pseudo")) {
  e <- learnErrors(files, pool=if(mode == "pooled") TRUE else "pseudo", multithread=FALSE, verbose=FALSE)
  write.table(matrix(sprintf("%.17g", e$err_out), 16), file.path(out, paste0(mode, "_err.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(e$trans, file.path(out, paste0(mode, "_trans.tsv")), sep="\t",
              row.names=FALSE, col.names=FALSE)
  writeLines(c(format(1e8, scientific=FALSE), length(e$err_in)), file.path(out, paste0(mode, "_info.txt")))
}
writeLines(c(paste("R", getRversion()), paste("dada2", packageVersion("dada2"))),
           file.path(out, "versions.txt"))
