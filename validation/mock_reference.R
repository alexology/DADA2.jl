suppressPackageStartupMessages(library(dada2))
root <- file.path("validation","generated","public")
f <- file.path(root,"r_filteredF.fastq.gz"); r <- file.path(root,"r_filteredR.fastq.gz")
stats <- filterAndTrim(file.path(root,"Mock_S280_L001_R1_001.fastq"),f,
                       file.path(root,"Mock_S280_L001_R2_001.fastq"),r,
                       truncLen=c(240,160),maxEE=c(2,2),maxN=0,truncQ=2,
                       rm.phix=FALSE,multithread=FALSE)
df <- derepFastq(f); dr <- derepFastq(r)
ef <- learnErrors(f,randomize=FALSE,multithread=FALSE,verbose=FALSE)
er <- learnErrors(r,randomize=FALSE,multithread=FALSE,verbose=FALSE)
fr <- dada(df,err=ef,multithread=FALSE,verbose=FALSE)
rr <- dada(dr,err=er,multithread=FALSE,verbose=FALSE)
merged <- mergePairs(fr,df,rr,dr,verbose=FALSE)
tab <- makeSequenceTable(list(merged))
tab <- tab[,nchar(colnames(tab)) %in% 250:256,drop=FALSE]
tab <- removeBimeraDenovo(tab,method="pooled",minFoldParentOverAbundance=1,verbose=FALSE)
write.table(data.frame(sequence=colnames(tab),count=as.integer(tab[1,])),file.path(root,"r_mock.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
write.table(getErrors(ef),file.path(root,"mock_errF.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
write.table(getErrors(er),file.path(root,"mock_errR.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
write.table(stats,file.path(root,"r_filter_stats.tsv"),sep="\t",quote=FALSE)
cat("Public mock: input=",stats[1,1]," filtered=",stats[1,2]," ASVs=",ncol(tab)," reads=",sum(tab),"\n")
