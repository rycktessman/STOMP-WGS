#SCRIPT TO GENERATE THE DISTANCE MATRIX

library(dplyr)
library(ggplot2)
library(tidyr)

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
