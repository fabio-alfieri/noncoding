# rm(ist=ls())
gc(full=T)

## produce SNV and CNA files from PCAWG ----
library(parallel)
library(plyr)
library(dplyr)

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

scratchMS <- "/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/"

columns <- read.table("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/patient_per_tumortype.tsv",
                      skip = 1)
tumor_types <- levels(factor(columns$V1))
tumor_types <- as.data.frame(cbind(do.call(rbind, str_split(tumor_types, "-")), tumor_types))

# for(tum in unique(tumor_types[duplicated(tumor_types$V1),]$V1)){
#   tumor_types[tumor_types$V1 == tum,]$tumor_types
#   snv_merged <- data.frame()
#   cna_merged <- data.frame()
#   for(i in tumor_types[tumor_types$V1 == tum,]$tumor_types){
#     snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
#     cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
#     snv_merged <- rbind.fill(snv_merged, snv)
#     cna_merged <- rbind.fill(cna_merged, cna)
#   }
#   write.table(cna_merged,
#               file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna/",tum,".merged_cna.tsv"),
#               sep = "\t")
#   write.table(snv_merged,
#               file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv/",tum,".merged_snv.tsv"),
#               sep = "\t")
# }

# exception Breast
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("BRCA-US", "BRCA-EU", "BRCA-UK")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Breast.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Breast.merged_snv.tsv"),
            sep = "\t")


# exception Colorectal
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("COAD-US", "READ-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Colorectal.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Colorectal.merged_snv.tsv"),
            sep = "\t")

# exception Liver
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("LIRI-JP", 
           "LIHC-US",
           "LINC-JP", "LICA-FR")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Liver.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Liver.merged_snv.tsv"),
            sep = "\t")

# exception Pancreas
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("PACA-CA", "PACA-AU", "PAEN-AU", "PAEN-IT"
           )){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Pancreas.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Pancreas.merged_snv.tsv"),
            sep = "\t")

# exception Brain
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("PBCA-DE", "GBM-US", "LGG-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Brain.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Brain.merged_snv.tsv"),
            sep = "\t")

# exception Kidney
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("RECA-EU"
           ,"KICH-US", "KIRC-US", "KIRP-US"
           )){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Kidney.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Kidney.merged_snv.tsv"),
            sep = "\t")

# exception Prostate
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("PRAD-CA", "EOPC-DE", "PRAD-UK", "PRAD-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Prostate.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Prostate.merged_snv.tsv"),
            sep = "\t")

# exception Ovary
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("OV-AU", "OV-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Ovary.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Ovary.merged_snv.tsv"),
            sep = "\t")

# exception Skin
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("MELA-AU", "SKCM-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Skin.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Skin.merged_snv.tsv"),
            sep = "\t")

# exception Lung
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("LUSC-US", "LUAD-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Lung.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Lung.merged_snv.tsv"),
            sep = "\t")


# exception Blood
snv_merged <- data.frame()
cna_merged <- data.frame()
for(i in c("MALY-DE", "CLLE-ES", "CMDI-UK", "LAML-KR", "DLBC-US")){
  snv <- read.table(file = paste0(scratchMS,"snv/",i,"_snv.tsv"))
  cna <- read.table(file = paste0(scratchMS,"cna/",i,"_cna.tsv"))
  snv_merged <- rbind.fill(snv_merged, snv)
  cna_merged <- rbind.fill(cna_merged, cna)
}
write.table(cna_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/cna.merged/Blood.merged_cna.tsv"),
            sep = "\t")
write.table(snv_merged,
            file = paste0("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/snv.merged/Blood.merged_snv.tsv"),
            sep = "\t")
