suppressPackageStartupMessages(library(dada2))
root <- file.path("test", "fixtures")
trans <- as.matrix(read.table(file.path(root,"error_transitions.tsv"),sep="\t"))
rownames(trans) <- paste0(rep(c("A","C","G","T"),each=4),"2",rep(c("A","C","G","T"),4))
colnames(trans) <- 0:(ncol(trans)-1)
direct <- matrix(0,16,ncol(trans))
for (source in 1:4) {
  rows <- (4*(source-1)+1):(4*source)
  total <- colSums(trans[rows,])
  for (destination in setdiff(1:4,source)) {
    row <- 4*(source-1)+destination
    df <- data.frame(q=0:(ncol(trans)-1),y=log10((trans[row,]+1)/total),total=total)
    df$y[!is.finite(df$y)] <- NA_real_
    model <- loess(y~q,df,weights=total,span=0.75,degree=2,
                   control=loess.control(surface="direct"))
    qrange <- range(df$q[is.finite(df$y)])
    prediction <- predict(model,data.frame(q=pmin(qrange[2],pmax(qrange[1],df$q))))
    direct[row,] <- pmax(1e-7,pmin(0.25,10^prediction))
  }
  direct[4*(source-1)+source,] <- 1-colSums(direct[setdiff(rows,4*(source-1)+source),])
}
write.table(direct,file.path(root,"error_loess_direct.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
interpolated <- loessErrfun(trans)
write.table(interpolated,file.path(root,"error_loess_r.tsv"),sep="\t",row.names=FALSE,col.names=FALSE)
writeLines(c(paste("R",getRversion()),paste("dada2",packageVersion("dada2")),
  "Quadratic weighted LOESS, span 0.75; direct surface and default interpolated surface."),
  file.path(root,"error_loess_versions.txt"))
