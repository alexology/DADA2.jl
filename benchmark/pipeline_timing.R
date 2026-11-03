# Paired-end pipeline in R dada2, timed by step, for comparison with
# benchmark/pipeline_timing.jl (same settings and outputs).
# Usage: Rscript benchmark/pipeline_timing.R R1.fastq.gz R2.fastq.gz out \
#          [trimLeftF trimLeftR truncLenF truncLenR maxEE threads]
suppressPackageStartupMessages(library(dada2))
args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 3) stop("Usage: pipeline_timing.R R1 R2 out [trimF trimR truncF truncR maxEE threads]")
num <- function(i, default) if(length(args) >= i) as.numeric(args[i]) else default
trim <- c(num(4, 0), num(5, 0)); trunc <- c(num(6, 0), num(7, 0)); maxee <- num(8, 2)
threads <- num(9, 4)
out <- file.path(args[3], "R"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
f <- file.path(out, "F.fastq.gz"); r <- file.path(out, "R.fastq.gz")
times <- list(); clock <- function() proc.time()[["elapsed"]]
t0 <- clock()
stats <- filterAndTrim(args[1], f, args[2], r, trimLeft=trim, truncLen=trunc,
                       maxEE=c(maxee, maxee), truncQ=2, maxN=0, multithread=threads)
times$filter <- clock() - t0; t0 <- clock()
ef <- learnErrors(f, multithread=threads); er <- learnErrors(r, multithread=threads)
times$learn <- clock() - t0; t0 <- clock()
df <- derepFastq(f); dr <- derepFastq(r)
times$derep <- clock() - t0; t0 <- clock()
af <- dada(df, err=ef, multithread=threads, verbose=FALSE)
ar <- dada(dr, err=er, multithread=threads, verbose=FALSE)
times$dada <- clock() - t0; t0 <- clock()
m <- mergePairs(af, df, ar, dr)
times$merge <- clock() - t0; t0 <- clock()
st <- makeSequenceTable(m)
nochim <- removeBimeraDenovo(st, method="consensus", multithread=threads)
times$chimeras <- clock() - t0
save_err <- function(e, name) write.table(matrix(sprintf("%.17g", getErrors(e)), 16), file.path(out, name),
                                          sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
save_err(ef, "err_F.tsv"); save_err(er, "err_R.tsv")
save_tab <- function(x, name) write.table(data.frame(sequence=names(x), abundance=unname(x)),
                                          file.path(out, name), sep="\t", row.names=FALSE, quote=FALSE)
save_tab(getUniques(af), "asv_F.tsv"); save_tab(getUniques(ar), "asv_R.tsv")
save_tab(getUniques(m), "merged.tsv")
save_tab(setNames(as.integer(st[1,]), colnames(st)), "table.tsv")
save_tab(setNames(as.integer(nochim[1,]), colnames(nochim)), "nochim.tsv")
writeLines(c(paste("reads_in", stats[1,1]), paste("reads_out", stats[1,2]),
             paste("rounds_F", length(ef$err_in)), paste("rounds_R", length(er$err_in)),
             sapply(names(times), function(k) paste0("time_", k, " ", round(times[[k]], 2)))),
           file.path(out, "summary.txt"))
cat(readLines(file.path(out, "summary.txt")), sep="\n")
