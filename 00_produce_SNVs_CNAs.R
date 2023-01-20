# rm(ist=ls())
gc(full=T)

## produce SNV and CNA files from PCAWG ----
library(parallel)
library(dplyr)
library(plyr)

pcawg.files <- list.files("~/Desktop_linux/mutation_compensation/data/clonal_analysis_PCAWG/")

cna <- mclapply(pcawg.files, mc.cores = 42, function(pt){
  fit <- readRDS(paste0("~/Desktop_linux/mutation_compensation/data/clonal_analysis_PCAWG/",pt,"/fit.rds"))
  cna <- fit$cna
  cna$ID <- pt
  return(cna)
})

cna_all <- do.call(rbind, cna)

snv <- mclapply(pcawg.files, mc.cores = 42, function(pt){
  fit <- readRDS(paste0("~/Desktop_linux/mutation_compensation/data/clonal_analysis_PCAWG/",pt,"/fit.rds"))
  snv <- fit$snv
  snv$ID <- pt
  return(snv)
})

columns <- data.frame()
for(i in 1:length(snv)){
  columns <- rbind(columns, cbind(i, dim(snv[[i]])[2], snv[[i]]$project_code[1]))
}

snv_final <- mclapply(levels(factor(columns$V3)), mc.cores = 4, function(tumor_type){
  do.call(rbind.fill, snv[columns$V3 == tumor_type])
})

names(snv_final) <- levels(factor(columns$V3))
# saveRDS(snv_final, file="snv_final.rsd")

mclapply(names(snv_final), mc.cores = 42, function(tumor_type){
  write.table(snv_final[[tumor_type]], 
              file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/",tumor_type,"_SNVs.tsv"),
              sep = "\t")
})


mclapply(names(snv_final), mc.cores = 42, function(tumor_type){
  write.table(cna_all[cna_all$ID %in% snv_final[[tumor_type]]$ID,],
              file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/",tumor_type,"_CNAs.tsv"),
              sep = "\t")
})
