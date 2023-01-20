# rm(list=ls())
gc(full=T)

suppressMessages({
  library(readxl)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(ggExtra)
  library(crayon)
  library(parallel)
})

scratchMS <- "/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/"
# columns <- read.table("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/patient_per_tumortype.tsv",
#                       skip = 1)

tumor_type <- "PACA-CA"

snv <- read.table(file = paste0(scratchMS,"snv/",tumor_type,"_snv.tsv"))
cna <- read.table(file = paste0(scratchMS,"cna/",tumor_type,"_cna.tsv"))

chr_info <-
  read.table("/home/ieo5099/Desktop_linux/mutation_compensation/data/misc/chr_info_h19.txt", header = TRUE)
chr_arms <-
  read.table(file = "/home/ieo5099/Desktop_linux/mutation_compensation/data/misc/cytoBand.txt", header = T)
chr_arms[, 2:3] <-
  apply(chr_arms[, 2:3] / fixed_bin_length, 2, as.integer)

for(chr in paste0("chr",1:22)){
  
}
