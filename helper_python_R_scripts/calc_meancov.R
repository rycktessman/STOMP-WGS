covg_files <- list.files(pattern = "\\.coverage")

summary_stats <- data.frame()
for(i in covg_files) {
  print(i)
  covg <- read.table(i)
  mean_depth <- mean(covg$V3)
  median_depth <- median(covg$V3)
  coverage10 <- sum(covg$V3>=10)/nrow(covg)
  coverage <- sum(covg$V3>0)/nrow(covg)
  
  info <- data.frame("id"=gsub("\\.coverage.*", "", i),
                     "mean_depth"=mean_depth,
                     "median_depth"=median_depth,
                     "prop_covg_over_10"=coverage10,
                     "prop_covg_any"=coverage)
  summary_stats <- rbind(summary_stats, info)
}

#ALSO GO THROUGH THOSE THAT FAILED QC
setwd("NoPassCov")
covg_files_failed <- list.files(pattern = "\\.coverage")
for(i in covg_files_failed) {
  print(i)
  covg <- read.table(i)
  mean_depth <- mean(covg$V3)
  median_depth <- median(covg$V3)
  coverage10 <- sum(covg$V3>=10)/nrow(covg)
  coverage <- sum(covg$V3>0)/nrow(covg)
  
  info <- data.frame("id"=gsub("\\.coverage.*", "", i),
                     "mean_depth"=mean_depth,
                     "median_depth"=median_depth,
                     "prop_covg_over_10"=coverage10,
                     "prop_covg_any"=coverage)
  summary_stats <- rbind(summary_stats, info)
}

write.csv(summary_stats, file="covg_summary.csv", row.names=F)
