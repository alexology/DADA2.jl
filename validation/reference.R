# Run from the repository root with Rscript validation/reference.R
suppressPackageStartupMessages(library(dada2))
out <- file.path("validation", "generated")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
root <- file.path("test", "fixtures")
f <- file.path(out, "forward.fastq.gz")
r <- file.path(out, "reverse.fastq.gz")
stats <- filterAndTrim(file.path(root,"sam1F.fastq.gz"),f,
                       file.path(root,"sam1R.fastq.gz"),r,
                       truncLen=c(240,160),maxEE=c(2,2),multithread=FALSE)
df <- derepFastq(f); dr <- derepFastq(r)
ef <- learnErrors(f,multithread=FALSE,randomize=FALSE,verbose=FALSE)
er <- learnErrors(r,multithread=FALSE,randomize=FALSE,verbose=FALSE)
fr <- dada(df,err=ef,multithread=FALSE,verbose=FALSE)
rr <- dada(dr,err=er,multithread=FALSE,verbose=FALSE)
write.table(getErrors(ef),file.path(out,"errors_forward.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
write.table(getErrors(er),file.path(out,"errors_reverse.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
write.table(data.frame(sequence=names(getUniques(fr)),count=unname(getUniques(fr))),
            file.path(out,"asvs_forward.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(sequence=names(getUniques(rr)),count=unname(getUniques(rr))),
            file.path(out,"asvs_reverse.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
merged <- mergePairs(fr,df,rr,dr,verbose=FALSE)
tab <- removeBimeraDenovo(makeSequenceTable(list(merged)),verbose=FALSE)
write.table(data.frame(sequence=colnames(tab),count=as.integer(tab[1,])),
            file.path(out,"asvs_final.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
cat("dada2 version:",as.character(packageVersion("dada2")),"\n")
print(stats)
cat("ASVs forward:",length(getUniques(fr)),"reverse:",length(getUniques(rr)),"\n")
