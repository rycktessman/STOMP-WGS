#SCRIPT TO GENERATE THE DISTANCE MATRIX

library(dplyr)
library(ggplot2)
library(tidyr)
#setwd("~/GitHub/STOMP_wgs/")

pair.dist <- function(id1, id2, data.dir, ext.filtered, ext.unfiltered) {
  select_cols = c("Chrom", "Position", "Ref", "Cons", "VarAllele",
                  "Reads1", "Reads2", "VarFreq", "Pvalue",
                  "Qual1", "Qual2", "Strands1", "Strands2")
  
  # rescue SNPs from unfiltered list, by two-way comparisons even if they didn't pass the depth and percent filters
  data1.unfiltered = read.table(file.path(data.dir, paste0(id1, ext.unfiltered)), header = T) %>% as_tibble() %>% dplyr::select(all_of(select_cols)) %>% tibble::add_column(type = "unfiltered")
  data1.filtered = read.table(file.path(data.dir, paste0(id1, ext.filtered)), header = T) %>% as_tibble() %>% dplyr::select(all_of(select_cols)) %>% tibble::add_column(type = "filtered")
  data2.unfiltered = read.table(file.path(data.dir, paste0(id2, ext.unfiltered)), header = T) %>% as_tibble() %>% dplyr::select(all_of(select_cols)) %>% tibble::add_column(type = "unfiltered")
  data2.filtered = read.table(file.path(data.dir, paste0(id2, ext.filtered)), header = T) %>% as_tibble() %>% dplyr::select(all_of(select_cols)) %>% tibble::add_column(type = "filtered")
  
  data1.exclude = low_reads_all[[id1]]
  data2.exclude = low_reads_all[[id2]]
  
  data1notin2.filtered =  data1.filtered %>% 
    dplyr::anti_join(data2.filtered, by = c("Position" = "Position"))
  # check unfiltered list from 2 to recover SNPs, and add them to data2 list
  data2.filtered.wrescued = data1notin2.filtered %>% 
    dplyr::inner_join(data2.unfiltered, by = c("Position" = "Position"), suffix = c(".data1", ".data2"), keep = T) %>%
    dplyr::select(ends_with(".data2")) %>%
    dplyr::rename_with(.fn = ~gsub(".data2", "", .), .cols = everything())
  if(nrow(data2.filtered) != 0) {
    data2.filtered.wrescued = data2.filtered.wrescued %>%
      dplyr::bind_rows(data2.filtered)
  }
  data2.filtered.wrescued = data2.filtered.wrescued %>%
    dplyr::rename_with(.fn = ~paste0(., ".", id2, sep = ""), .cols = -c(Chrom, Position, Ref))
  
  
  data2notin1.filtered =  data2.filtered %>% 
    dplyr::anti_join(data1.filtered, by = c("Position" = "Position"))
  # check unfiltered list from 2 to recover SNPs, and add them to data2 list
  data1.filtered.wrescued = data2notin1.filtered %>% 
    dplyr::inner_join(data1.unfiltered, by = c("Position" = "Position"), suffix = c(".data2", ".data1"), keep = T) %>%
    dplyr::select(ends_with(".data1")) %>%
    dplyr::rename_with(.fn = ~gsub(".data1", "", .), .cols = everything())
  if(nrow(data1.filtered) != 0) {
    data1.filtered.wrescued = data1.filtered.wrescued %>%
      dplyr::bind_rows(data1.filtered)
  }
  data1.filtered.wrescued = data1.filtered.wrescued %>%
    dplyr::rename_with(.fn = ~paste0(., ".", id1, sep = ""), .cols = -c(Chrom, Position, Ref))
  
  data12 = data1.filtered.wrescued %>% 
    dplyr::full_join(data2.filtered.wrescued, by = c("Position" = "Position")) %>%
    dplyr::select(Chrom.x, Position, Ref.x, 
                  paste0("Cons.", id1),
                  paste0("VarAllele.", id1),
                  paste0("Reads1.", id1),
                  paste0("Reads2.", id1),
                  paste0("VarFreq.", id1),
                  paste0("type.", id1),
                  Chrom.y, Ref.y, 
                  paste0("Cons.", id2),
                  paste0("VarAllele.", id2),
                  paste0("Reads1.", id2),
                  paste0("Reads2.", id2),
                  paste0("VarFreq.", id2),
                  paste0("type.", id2)) %>%
    dplyr::mutate(Chrom.x = case_when(is.na(Chrom.x) ~ "MTB_anc",
                                      TRUE ~ Chrom.x)) %>%
    dplyr::mutate(Ref.x = case_when(is.na(Ref.x) ~ Ref.y,
                                    TRUE ~ Ref.x)) %>%
    dplyr::mutate_at(vars(matches("VarFreq")), ~ as.numeric(sub("%", "", .))) %>%
    dplyr::rename(Chrom = Chrom.x, Ref = Ref.x) %>%
    dplyr::select(-Chrom.y, -Ref.y) %>%
    dplyr::arrange(Position) %>%
    dplyr::mutate(comparison = case_when(!is.na(get(paste0("VarAllele.", id1))) & is.na(get(paste0("VarAllele.", id2))) ~ "unique",
                                         is.na(get(paste0("VarAllele.", id1))) & !is.na(get(paste0("VarAllele.", id2))) ~ "unique",
                                         ((!is.na(get(paste0("VarAllele.", id1)))) & 
                                            (!is.na(get(paste0("VarAllele.", id2)))) &
                                            (get(paste0("VarAllele.", id1)) == get(paste0("VarAllele.", id2)))) ~ "common",
                                         ((!is.na(get(paste0("VarAllele.", id1)))) & 
                                            (!is.na(get(paste0("VarAllele.", id2)))) &
                                            (get(paste0("VarAllele.", id1)) != get(paste0("VarAllele.", id2)))) ~ "unique",
                                         TRUE ~ "other")) %>%
    dplyr::select(Chrom, Position, Ref, comparison, everything())
  
  data1.pos.rescued = data12 %>% dplyr::filter(get(paste0("type.", id1)) %in% c("filtered", "unfiltered")) %>% dplyr::pull(Position)
  data2.pos.rescued = data12 %>% dplyr::filter(get(paste0("type.", id2)) %in% c("filtered", "unfiltered")) %>% dplyr::pull(Position)
  data1.nrescued = data1.pos.rescued %>% length()
  data2.nrescued = data2.pos.rescued %>% length()
  intersect.nrescued = data12 %>% dplyr::filter(comparison == "common") %>% dplyr::pull(Position) %>% length()
  #count SNPs if position gets called as a SNP in both sequences, but most frequent non-reference allele differs
  diff.nrescued.orig = data12 %>% 
    dplyr::filter(comparison == "unique" |
                    get(paste0("VarAllele.", id1)) != get(paste0("VarAllele.", id2))) %>% 
    dplyr::pull(Position) %>% length()
  diff.nrescued.exclusions = data12 %>% 
    dplyr::filter(comparison == "unique" |
                    get(paste0("VarAllele.", id1)) != get(paste0("VarAllele.", id2))) %>% 
    dplyr::filter(!(Position %in% data1.exclude) & !(Position %in% data2.exclude)) %>% 
    dplyr::pull(Position) %>% length()
  
  result = data.frame(id1 = id1, 
                      id2 = id2, 
                      id1.nrescued = data1.nrescued,
                      id2.nrescued = data2.nrescued,
                      intersect.nrescued = intersect.nrescued,
                      diff.nrescued.orig = diff.nrescued.orig,
                      diff.nrescued = diff.nrescued.exclusions)
  result
}

check_symetry <- function(matrix) {
  for(i in 1:nrow(matrix)) {
    for(j in (i+1):ncol(matrix)) {
      a = matrix[i, j] - matrix[j, i]
      cat(i, j, matrix[i, j], matrix[j, i],"\n")
      #if(a != 0) {cat(i, j, "\n")}
    }
  }
}

data.dir = "fq_files/"
path_out = "fq_files/"
ids.wSNPs = sub(".fSNPs.annoF.delF.densF", "", list.files(data.dir, pattern=".fSNPs.annoF.delF.densF"))
#remaining IDs failed QC

#add list of low-read depth positions to be excluded from a given isolate's SNPs
load(paste0(data.dir, "low_read_positions.Rda"))

ids.pairs = expand.grid(ids.wSNPs, ids.wSNPs, stringsAsFactors = F) %>% as_tibble() %>% dplyr::filter(Var1 != Var2)

ext.filtered <- ".fSNPs.annoF.delF.densF" #with annotations filtering, indel filtering, window-based filtering
ext.unfiltered <- ".snp" #for rescue SNPs - without any filtering

#server version
if(FALSE) {
  distances = parallel::mclapply(1:nrow(ids.pairs),
                                 function(r) {
                                   result = pair.dist(ids.pairs[r, 1], ids.pairs[r, 2], 
                                                      data.dir = data.dir,
                                                      ext.filtered = ext.filtered,
                                                      ext.unfiltered = ext.unfiltered)
                                   result
                                 },
                                 mc.cores = 40) 
  distances_rbind = lapply(distances,"[[",1)
  distances_rbind = do.call("rbind", distances_rbind)
}

#loop version
distances_rbind <- data.frame()
for(r in 1:nrow(ids.pairs)) {
  print(r)
  result <- pair.dist(ids.pairs[[r, 1]], ids.pairs[[r, 2]],
                      data.dir=data.dir,
                      ext.filtered=ext.filtered,
                      ext.unfiltered=ext.unfiltered)
  distances_rbind <- bind_rows(distances_rbind, result)
}

#excluding low read depths
distances_rbind.diff.nrescued = distances_rbind %>% as_tibble() %>% 
  dplyr::rename(sample1 = id1, sample2 = id2, dist = diff.nrescued) %>% 
  dplyr::filter(sample1 != sample2) %>%
  dplyr::arrange(dist) %>%
  dplyr::select(sample1, dist, sample2)

#convert to matrix
dist_append <- data.frame("sample1"=sort(unique(distances_rbind.diff.nrescued$sample1)), 
                          "dist"=0, 
                          "sample2"=sort(unique(distances_rbind.diff.nrescued$sample2)))
distances_rbind.diff.nrescued.matrix = rbind(distances_rbind.diff.nrescued,
                                             dist_append) %>%
  tidyr::pivot_wider(names_from = sample1, values_from = dist) %>%
  tibble::column_to_rownames(var = "sample2")
distances_rbind.diff.nrescued.matrix <- distances_rbind.diff.nrescued.matrix[order(row.names(distances_rbind.diff.nrescued.matrix)) ,]
distances_rbind.diff.nrescued.matrix <- distances_rbind.diff.nrescued.matrix[, order(colnames(distances_rbind.diff.nrescued.matrix))]
print(isSymmetric(as.matrix(distances_norescue_rbind.diff.wlowreads.matrix)))


write.csv(distances_rbind.diff.nrescued, file = paste0(path_out, "pairs-dist-densityfiltered_wrescueSNPs.csv"), quote = F, row.names = FALSE)
write.csv(distances_rbind.diff.nrescued.matrix, file = paste0(path_out, "distmatrix-densityfiltered_wrescueSNPs.csv"), quote = F, row.names = T)
write.table(distances_rbind.diff.nrescued, file = paste0(path_out, "pairs-dist-densityfiltered_wrescueSNPs.tsv"), quote = F, row.names = FALSE, sep="\t")
write.table(distances_rbind.diff.nrescued.matrix, file = paste0(path_out, "distmatrix-densityfiltered_wrescueSNPs.tsv"), quote = F, row.names = T, col.names=NA, sep="\t")

#including low read depths 
distances_rbind.diff.wlowreads = distances_rbind %>% as_tibble() %>% 
  dplyr::rename(sample1 = id1, sample2 = id2, dist = diff.nrescued.orig) %>% 
  dplyr::filter(sample1 != sample2) %>%
  dplyr::arrange(dist) %>%
  dplyr::select(sample1, dist, sample2)

#convert to matrix
dist_append <- data.frame("sample1"=sort(unique(distances_rbind.diff.wlowreads$sample1)), 
                          "dist"=0, 
                          "sample2"=sort(unique(distances_rbind.diff.wlowreads$sample2)))
distances_rbind.diff.wlowreads.matrix = rbind(distances_rbind.diff.wlowreads,
                                             dist_append) %>%
  tidyr::pivot_wider(names_from = sample1, values_from = dist) %>%
  tibble::column_to_rownames(var = "sample2")
distances_rbind.diff.wlowreads.matrix <- distances_rbind.diff.wlowreads.matrix[order(row.names(distances_rbind.diff.wlowreads.matrix)) ,]
distances_rbind.diff.wlowreads.matrix <- distances_rbind.diff.wlowreads.matrix[, order(colnames(distances_rbind.diff.wlowreads.matrix))]
print(isSymmetric(as.matrix(distances_norescue_rbind.diff.wlowreads.matrix)))

write.csv(distances_rbind.diff.wlowreads, file = paste0(path_out, "pairs-dist-densityfiltered_wrescueSNPs_wlowreads.csv"), quote = F, row.names = FALSE)
write.csv(distances_rbind.diff.wlowreads.matrix, file = paste0(path_out, "distmatrix-densityfiltered_wrescueSNPs_wlowreads.csv"), quote = F, row.names = T)
write.table(distances_rbind.diff.wlowreads, file = paste0(path_out, "pairs-dist-densityfiltered_wrescueSNPs_wlowreads.tsv"), quote = F, row.names = FALSE, sep="\t")
write.table(distances_rbind.diff.wlowreads.matrix, file = paste0(path_out, "distmatrix-densityfiltered_wrescueSNPs_wlowreads.tsv"), quote = F, row.names = T, col.names=NA, sep="\t")


pairs = list()
pairs[[1]] = c("1A", "56124_S488_L005", "STP_56124A_S217_L004")
pairs[[2]] = c("1A", "56124_S488_L005", "STP_56124B_S236_L004")
pairs[[3]] = c("1A", "STP_56124A_S217_L004", "STP_56124B_S236_L004")
pairs[[4]] = c("2A", "57875B_S485_L005", "57875A_S478_L005")
pairs[[5]] = c("3A", "STP61302_S637_L006", "STP61313_S627_L006")
pairs[[6]] = c("4A", "STP62336_S655_L006", "STP62359_S649_L006")
pairs[[7]] = c("5A", "STP62440_S673_L006", "STP62498_S642_L006")
pairs[[8]] = c("6A", "5872_S85", "STP_35872_S261_L004")
pairs[[9]] = c("7A", "5241_S96", "STP_35241_S213_L004")
pairs[[10]] = c("8A", "37316_S513_L005", "STP_37316_S222_L004")
pairs[[11]] = c("1B", "60261_S595_L006", "STP60261-S29")
pairs[[12]] = c("2B", "STP63747-S53", "STP63409_S658_L006")
pairs[[13]] = c("3B", "STP62439_S665_L006", "STP62499-S51")
pairs[[14]] = c("4B", "STP69976-S44", "STP70469-S31")
pairs[[15]] = c("5B", "STP68461-S41", "STP68351-S46")
pairs[[16]] = c("6B", "57875B_S485_L005", "STP57875-S4")
pairs[[17]] = c("6B", "57875A_S478_L005", "STP57875-S4")
pairs[[18]] = c("7B", "STP63533-S49", "STP63533_S660_L006")
pairs[[19]] = c("1C", "STP61695-S11", "FT100020573_L01_STP_61695")
pairs[[20]] = c("2C","STP_56121_S238_L004", "FT100020573_L01_STP_56121")
pairs[[21]] = c("3C","STP_37679_S207_L004", "FT100020573_L01_STP_37679")
pairs[[22]] = c("4C","STP62498_S642_L006", "FT100020573_L01_STP_62498")
pairs[[23]] = c("5C","STP63114_S648_L006", "FT100020573_L01_STP_63114")
pairs[[24]] = c("6C","4706_S101", "FT100020573_L01_STP_34706")

for(i in 1:length(pairs)) {
  temp = distances_rbind.diff.nrescued %>%
    dplyr::filter(sample1 == pairs[[i]][2] & sample2 == pairs[[i]][3]) %>%
    as.data.frame()
  cat(pairs[[i]][1], "\t", temp[1, "sample1"], "\t", temp[1, "dist"], "\t", temp[1, "sample2"], "\n")
}



