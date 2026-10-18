data <- read.delim(file.path("validation","error_learning_results.tsv"))
stopifnot(nrow(data)==210, length(unique(data$method))==7, all(table(data$scenario,data$method)==5))
metrics <- c("substitution_l1","log10_rmse","predicted_error_rate","true_error_rate",
             "error_ratio","precision","recall_ge2","abundance_tv")
summary <- aggregate(data[metrics],data[c("scenario","method")],
                     function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE))
write.table(summary,file.path("validation","error_learning_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
directory <- file.path("docs","figures")
dir.create(directory,showWarnings=FALSE)
draw <- function() {
  par(mfrow=c(1,2),mar=c(7,4.5,3,1))
  q <- 25:40
  oracle <- as.matrix(read.table(file.path("validation","generated","science","calibrated_11","oracle.tsv"),sep="\t"))
  methods <- c("julia_linear","julia_loess","julia_loess_raw","julia_loess_raw_shared")
  models <- lapply(methods,function(m) as.matrix(read.table(file.path("validation","generated","error_learning",paste0("calibrated_11_",m,".tsv")),sep="\t")))
  colors <- c("#777777","#d95f02","#0072b2","#009e73")
  y <- cbind(oracle[3,q+1],sapply(models,function(m) m[3,q+1]))
  matplot(q,y,type="l",log="y",lty=c(1,2,2,1,1),lwd=2,
          col=c("black",colors),xlab="Quality score",ylab="A to G substitution probability",
          main="Known model vs estimated rates\nCalibrated simulation, seed 11")
  legend("bottomleft",c("Truth","Linear, mean Q","LOESS, mean Q","LOESS, raw Q","LOESS, raw Q, shared"),
         col=c("black",colors),lty=c(1,2,2,1,1),lwd=2,cex=.72,bty="n")
  scenarios <- c("calibrated","rare_snv","miscalibrated","indels","chimeras")
  method <- c("julia_linear","julia_loess","julia_loess_raw","julia_loess_raw_shared")
  ratios <- sapply(scenarios,function(s) sapply(method,function(m) mean(data$error_ratio[data$scenario==s & data$method==m])))
  barplot(ratios,beside=TRUE,log="y",col=colors,names.arg=c("Calibrated","Rare SNV","5x errors","Indels","Chimeras"),
          las=2,ylab="Estimated / true error rate (log scale)",main="Error calibration\nMean across five seeds")
  abline(h=1,lty=3)
}
svg(file.path(directory,"error-learning.svg"),width=11,height=6)
draw(); dev.off()
png(file.path(directory,"error-learning.png"),width=2200,height=1200,res=200)
draw(); dev.off()
