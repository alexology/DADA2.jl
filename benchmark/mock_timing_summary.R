paths <- c("mock_timing_R_1.tsv","mock_timing_Julia_1.tsv","mock_timing_R_4.tsv","mock_timing_Julia_4.tsv")
data <- do.call(rbind,lapply(paths,function(path) read.delim(file.path("benchmark",path))))
stopifnot(nrow(data)==60,all(data$reads_in==4779),all(data$filtered==4314),all(data$asvs==20),all(data$final_reads==4269))
reference <- read.delim(file.path("validation","generated","public","r_mock.tsv"))
reference <- reference[order(reference$sequence),]
data$identical_counts_to_R <- FALSE
data$count_l1_difference <- NA_integer_
for(i in seq_len(nrow(data))) {
  row <- data[i,]
  prefix <- if(row$implementation=="R") "r" else "julia"
  path <- file.path("validation","generated","mock_timing",paste0(prefix,"_",row$threads,"_",row$mode,"_",row$repetition,".tsv"))
  actual <- read.delim(path); actual <- actual[order(actual$sequence),]
  same_sequences <- identical(actual$sequence,reference$sequence)
  data$identical_counts_to_R[i] <- same_sequences && all(actual$count==reference$count)
  if(same_sequences) data$count_l1_difference[i] <- sum(abs(actual$count-reference$count))
}
stopifnot(all(data$identical_counts_to_R[data$mode=="shared_errors"]))
write.table(data,file.path("benchmark","mock_timing_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
warm <- data[data$phase=="warm",]
summary <- do.call(rbind,lapply(split(warm,interaction(warm$implementation,warm$threads,warm$mode,drop=TRUE)),function(group) {
  stopifnot(nrow(group)==5)
  result <- group[1,c("implementation","threads","mode")]
  result$repetitions <- nrow(group)
  for(name in c("filter_seconds","derep_seconds","learning_seconds","inference_seconds","merge_seconds","community_seconds","total_seconds"))
    result[[paste0("median_",name)]] <- median(group[[name]])
  result$min_total_seconds <- min(group$total_seconds); result$max_total_seconds <- max(group$total_seconds)
  result$identical_counts_to_R <- all(group$identical_counts_to_R)
  result
}))
rownames(summary) <- NULL
summary <- summary[order(summary$threads,summary$mode,summary$implementation),]
write.table(summary,file.path("benchmark","mock_timing_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(summary[c("implementation","threads","mode","median_total_seconds","min_total_seconds","max_total_seconds","identical_counts_to_R")],row.names=FALSE)
cat("\nFirst result after entering script (interpreter startup excluded):\n")
print(data[!is.na(data$first_script_seconds),c("implementation","threads","first_script_seconds")],row.names=FALSE)
