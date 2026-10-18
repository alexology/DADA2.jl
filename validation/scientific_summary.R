root <- "validation"
data <- read.delim(file.path(root,"scientific_results.tsv"))
stopifnot(nrow(data)==150, all(table(data$case)==5))
scenario_order <- c("no_errors","calibrated","rare_snv","miscalibrated","indels","chimeras")
method_order <- c("julia_oracle","julia_learned","julia_r_errors","r_oracle","r_learned")
metrics <- c("precision","recall","recall_ge2","abundance_tv","spurious_read_fraction")
summary <- data.frame()
for (scenario in scenario_order) for (method in method_order) {
  rows <- data[data$scenario==scenario & data$method==method,]
  stopifnot(nrow(rows)==5)
  entry <- data.frame(scenario=scenario,method=method,n=nrow(rows))
  for (metric in metrics) {
    entry[[paste0(metric,"_mean")]] <- mean(rows[[metric]])
    entry[[paste0(metric,"_min")]] <- min(rows[[metric]])
    entry[[paste0(metric,"_max")]] <- max(rows[[metric]])
  }
  summary <- rbind(summary,entry)
}
write.table(summary,file.path(root,"scientific_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
dir.create("docs/figures",recursive=TRUE,showWarnings=FALSE)
draw <- function() {
  par(mfrow=c(3,1),mar=c(3,4.5,2,1),oma=c(3,0,1,0),las=1)
  shown <- c("julia_oracle","julia_learned","r_oracle","r_learned")
  colors <- c("#0072B2","#CC79A7","#009E73","#D55E00")
  labels <- c("No errors","Calibrated","Rare SNV","5x error","Indels","Chimeras")
  for (metric in c("precision","recall_ge2","abundance_tv")) {
    scale <- if(metric=="abundance_tv") 1 else 100
    ylim <- if(metric=="abundance_tv") c(0,0.06) else if(metric=="precision") c(70,102) else c(50,102)
    title <- switch(metric,precision="Exact-sequence precision (%)",recall_ge2="Sensitivity excluding biological singletons (%)",
                    abundance_tv="Relative-abundance total variation distance")
    plot(NA,xlim=c(0.5,6.5),ylim=ylim,xaxt="n",xlab="",ylab=title,bty="l")
    axis(1,at=1:6,labels=labels,cex.axis=.9)
    abline(h=pretty(ylim),col="#EEEEEE",lty=1)
    for (j in seq_along(shown)) {
      rows <- summary[summary$method==shown[j],]
      rows <- rows[match(scenario_order,rows$scenario),]
      x <- (1:6)+(j-2.5)*0.08
      y <- rows[[paste0(metric,"_mean")]]*scale
      lo <- rows[[paste0(metric,"_min")]]*scale
      hi <- rows[[paste0(metric,"_max")]]*scale
      lines(x,y,col=colors[j],lty=j,lwd=1.5)
      points(x,y,col=colors[j],pch=15+j,cex=1.1)
      varying <- hi-lo > 1e-9
      if(any(varying)) arrows(x[varying],lo[varying],x[varying],hi[varying],angle=90,code=3,length=.04,col=colors[j])
    }
    if(metric=="precision") legend("bottomleft",c("Julia known errors","Julia learned errors","R known errors","R learned errors"),
                 col=colors,lty=1:4,pch=16:19,ncol=2,bty="n",cex=.85)
  }
  mtext("Five seeded communities per condition; markers = means, bars = observed ranges (not confidence intervals)",
        outer=TRUE,side=1,line=1,cex=.8)
}
svg("docs/figures/scientific-validation.svg",width=10,height=10); draw(); dev.off()
png("docs/figures/scientific-validation.png",width=1800,height=1800,res=180); draw(); dev.off()
cat("Validated",nrow(data),"method/sample measurements; summary and scientific figures written.\n")
