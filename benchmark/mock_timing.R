script_start <- proc.time()[["elapsed"]]
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
threads <- as.integer(args[[1]])
Sys.setenv(OMP_NUM_THREADS=threads,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1)
suppressPackageStartupMessages(library(dada2))
load_seconds <- proc.time()[["elapsed"]]-script_start
root <- file.path("validation","generated","public")
countsdir <- file.path("validation","generated","mock_timing")
dir.create(countsdir,recursive=TRUE,showWarnings=FALSE)
shared <- lapply(c("mock_errF.tsv","mock_errR.tsv"),function(name)
  as.matrix(read.table(file.path(root,name),sep="\t")))
timed <- function(f) {
  start <- proc.time()[["elapsed"]]
  value <- f()
  list(value=value,seconds=proc.time()[["elapsed"]]-start)
}
pipeline <- function(dir,mode) {
  f <- file.path(dir,"F.fastq.gz"); r <- file.path(dir,"R.fastq.gz")
  filtering <- timed(function() filterAndTrim(file.path(root,"Mock_S280_L001_R1_001.fastq"),f,
    file.path(root,"Mock_S280_L001_R2_001.fastq"),r,truncLen=c(240,160),maxEE=c(2,2),
    maxN=0,truncQ=2,minLen=20,rm.phix=FALSE,multithread=FALSE,OMP=threads>1,
    compress=TRUE,verbose=FALSE))
  dereplication <- timed(function() list(derepFastq(f,verbose=FALSE),derepFastq(r,verbose=FALSE)))
  df <- dereplication$value[[1]]; dr <- dereplication$value[[2]]
  learning_seconds <- 0
  iterations_F <- iterations_R <- 0
  if(mode=="shared_errors") {
    ef <- shared[[1]]; er <- shared[[2]]
  } else {
    estimate <- function(d) {
      iterations <- 0L
      estimator <- function(trans) { iterations <<- iterations+1L; loessErrfun(trans) }
      err <- learnErrors(d,nbases=1e8,randomize=FALSE,multithread=if(threads==1) FALSE else threads,
        MAX_CONSIST=10,OMEGA_C=0,errorEstimationFunction=estimator,verbose=0)
      list(errors=getErrors(err),iterations=iterations)
    }
    learning <- timed(function() list(estimate(df),estimate(dr)))
    ef <- learning$value[[1]]$errors; er <- learning$value[[2]]$errors
    iterations_F <- learning$value[[1]]$iterations; iterations_R <- learning$value[[2]]$iterations
    learning_seconds <- learning$seconds
  }
  inference <- timed(function() list(
    dada(df,err=ef,multithread=if(threads==1) FALSE else threads,verbose=FALSE),
    dada(dr,err=er,multithread=if(threads==1) FALSE else threads,verbose=FALSE)))
  fr <- inference$value[[1]]; rr <- inference$value[[2]]
  merging <- timed(function() mergePairs(fr,df,rr,dr,minOverlap=12,maxMismatch=0,verbose=FALSE))
  community <- timed(function() {
    tab <- makeSequenceTable(list(merging$value))
    tab <- tab[,nchar(colnames(tab)) %in% 250:256,drop=FALSE]
    removeBimeraDenovo(tab,method="pooled",minFoldParentOverAbundance=1,verbose=FALSE)
  })
  list(table=community$value,stats=filtering$value,iterations=c(iterations_F,iterations_R),
    times=c(filtering$seconds,dereplication$seconds,learning_seconds,inference$seconds,merging$seconds,community$seconds))
}
header <- paste(c("implementation","threads","mode","repetition","phase","filter_seconds","derep_seconds",
  "learning_seconds","inference_seconds","merge_seconds","community_seconds","total_seconds",
  "script_load_seconds","first_script_seconds","reads_in","filtered","asvs","final_reads",
  "error_iterations_F","error_iterations_R","converged"),collapse="\t")
writeLines(header,args[[2]])
for(mode in c("full_mean","shared_errors")) for(rep in 0:5) {
  gc()
  dir <- tempfile("dada2-timing-"); dir.create(dir)
  result <- tryCatch(timed(function() pipeline(dir,mode)),finally=unlink(dir,recursive=TRUE))
  tab <- result$value$table; stats <- result$value$stats
  first_seconds <- if(mode=="full_mean" && rep==0) proc.time()[["elapsed"]]-script_start else "NA"
  stopifnot(stats[1,1]==4779,stats[1,2]==4314,ncol(tab)==20,sum(tab)==4269)
  values <- c("R",threads,mode,rep,if(rep==0) "first_call" else "warm",result$value$times,
    result$seconds,load_seconds,first_seconds,stats[1,1],stats[1,2],ncol(tab),sum(tab),result$value$iterations,"NA")
  cat(paste(values,collapse="\t"),"\n",file=args[[2]],append=TRUE,sep="")
  write.table(data.frame(sequence=colnames(tab),count=as.integer(tab[1,])),
    file.path(countsdir,paste0("r_",threads,"_",mode,"_",rep,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  cat("R",threads,"threads",mode,"rep",rep,":",round(result$seconds,3),"s\n"); flush.console()
}
writeLines(c(paste("R",getRversion()),paste("dada2",packageVersion("dada2"))),
  file.path("benchmark","mock_timing_R_versions.txt"))
writeLines(capture.output(print(warnings())),file.path("benchmark",paste0("mock_timing_R_",threads,"_warnings.txt")))
