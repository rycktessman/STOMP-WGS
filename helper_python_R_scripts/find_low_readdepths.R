#using coverage files, find positions with low read depths
#generate these in a separate RDA file since the covg files are large and take a long time to load from the server
#exclude these positions (from a given strain only) in the distance matrix

library(data.table)

threshold <- 3 #tag positions with < 3 reads - because with 3 reads or more they'll show up in rescue SNPs

data.dir = "fq_files/"

ids.wSNPs = sub(".fSNPs.annoF.delF.densF", "", list.files(data.dir, pattern=".densF"))

low_reads_all <- list()
for(i in ids.wSNPs[1:length(ids.wSNPs)]) {
  print(i)
  mpileup <- fread(paste0(data.dir, i, ".mpileup"))
  low_reads <- mpileup[V4<threshold, "V2"][[1]]
  no_reads <-  (1:4411532)[!(1:4411532 %in% mpileup$V2)] #if 0 reads, doesn't show up in mpileup file
  low_reads_all[[i]] <- c(low_reads, no_reads)
}

save(low_reads_all, file=paste0(data.dir, "low_read_positions.Rda"))
