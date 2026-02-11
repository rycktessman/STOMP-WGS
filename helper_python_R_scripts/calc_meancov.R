#path <- "fq_files_TRAC_test/"
#path <- Sys.getenv('path')
#setwd(path)
covg_files <- list.files(pattern = "\\.coverage")

#CONFIRMED THAT THIS VERSION WORKS - DON'T NEED SLOWER ALT VERSION BELOW
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

#slower version that produces the same thing (from Inaki's pipeline)
if(FALSE) {
  for(i in covg_files) {
    print(i)
    covg <- read.table(i)
    positions <- 0
    coverage <- 0
    covered_positions <- 0
    cov_list <- c()
    for(j in 1:nrow(covg)) {
      ref <- covg[j, 1]
      pos <- covg[j, 2]
      cov <- covg[j, 3]
      cov <- as.integer(cov)
      positions <- positions + 1
      coverage <- coverage + cov
      if(cov >= 10) {
        covered_positions <- covered_positions + 1
      }
      cov_list <- c(cov_list, cov)
    }
    mean_cov <- coverage/positions
    genome_coverage <- covered_positions/positions
    
  }
}


write.csv(summary_stats, file="covg_summary.csv", row.names=F)
