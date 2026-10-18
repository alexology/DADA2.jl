suppressPackageStartupMessages(library(dada2))
root <- file.path("validation","generated","science")
manifest <- read.delim(file.path(root,"manifest.tsv"))
write_result <- function(result,path) {
  tab <- makeSequenceTable(list(result))
  tab <- removeBimeraDenovo(tab,method="pooled",minFoldParentOverAbundance=1,verbose=FALSE)
  write.table(data.frame(sequence=colnames(tab),count=as.integer(tab[1,])),path,
              sep="\t",quote=FALSE,row.names=FALSE)
}
status <- data.frame(case=character(),mode=character(),error=character())
for (name in manifest$case) {
  dir <- file.path(root,name)
  cat("R reference:",name,"\n"); flush.console()
  d <- derepFastq(file.path(dir,"reads.fastq"),verbose=FALSE)
  oracle <- as.matrix(read.table(file.path(dir,"oracle.tsv"),sep="\t"))
  for (mode in c("oracle","learned")) {
    problem <- tryCatch({
      err <- if(mode=="oracle") oracle else getErrors(learnErrors(
        file.path(dir,"reads.fastq"),nbases=1e8,randomize=FALSE,multithread=FALSE,verbose=FALSE))
      result <- dada(d,err=err,multithread=FALSE,verbose=FALSE)
      write_result(result,file.path(dir,paste0("r_",mode,".tsv")))
      if(mode=="learned") write.table(err,file.path(dir,"r_errors.tsv"),sep="\t",col.names=FALSE,row.names=FALSE)
      ""
    },error=function(e) conditionMessage(e))
    status <- rbind(status,data.frame(case=name,mode=mode,error=problem))
  }
}
write.table(status,file.path(root,"r_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(paste("R",getRversion()),paste("dada2",packageVersion("dada2"))),file.path(root,"r_versions.txt"))
