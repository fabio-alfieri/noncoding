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

setwd("~/mountHD/noncoding/")

scratchMS <- "/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/"
fixed_bin_length <- 1000000

results_table_path <- "results/01_binLevel_noIntron/"

columns <- read.table("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/patient_per_tumortype.tsv",
                      skip = 1)
if(F){
  tumor_types <- levels(factor(columns$V1))
  tumor_types <- as.data.frame(cbind(do.call(rbind, str_split(tumor_types, "-")), tumor_types))
  tumor_types <- tumor_types$tumor_types
}

# tumor_types <- c("BLCA-US","BOCA-UK","BRCA.merged","BTCA-SG","CESC-US","CLLE-ES","CMDI-UK",
# "COADREAD.merged","DLBC-US","EOPC-DE","ESAD-UK","GACA-CN","GBM-US","HNSC-US","KICH-US","KIRC-US",
# "KIRP-US","LAML-KR","LGG-US","LICA-FR","LIHC-US","LINC-JP","LIRI-JP","LUAD-US","LUSC-US",
# "MALY-DE","MELA-AU","ORCA-IN","OV.merged","PACA.merged","PAEN.merged",
# "PBCA-DE","PRAD-CA","PRAD.merged","RECA-EU","SARC-US","SKCM-US","STAD-US",
# "THCA-US","UCEC-US")
# 
tumor_types <- paste0(c(
                 # "Breast",
                 "Colorectal",
                 "Brain",
                 "Kidney",
                 "Liver",
                 "Pancreas",
                 "Lung",
                 # "Prostate",
                 "Ovary",
                 "Skin",
                 "Blood"
                 ), ".merged")
merged <- T

mclapply(tumor_types, mc.cores = 9, function(tumor_type){
# for(tumor_type in tumor_types){
  print(tumor_type)
  if(merged){
    snv <- read.table(file = paste0(scratchMS,"snv.merged/",tumor_type,"_snv.tsv"))
    cna <- read.table(file = paste0(scratchMS,"cna.merged/",tumor_type,"_cna.tsv"))
  }else{
    snv <- read.table(file = paste0(scratchMS,"snv/",tumor_type,"_snv.tsv"))
    cna <- read.table(file = paste0(scratchMS,"cna/",tumor_type,"_cna.tsv"))
  }
  
  common_patients <- levels(factor(snv$ID))[levels(factor(snv$ID)) %in% levels(factor(cna$ID))]
  
  chr_info <-
    read.table("/home/ieo5099/Desktop_linux/mutation_compensation/data/misc/chr_info_h19.txt", header = TRUE)
  chr_arms <-
    read.table(file = "/home/ieo5099/Desktop_linux/mutation_compensation/data/misc/cytoBand.txt", header = T)
  chr_arms[, 2:3] <-
    apply(chr_arms[, 2:3] / 1000000, 2, as.integer)
  
  # Start loop for cancer type ----
  # for(chr in paste0("chr",1:22)){
  mclapply(paste0("chr", 1:22), mc.cores = 2, function(chr){
    
    # filter for chromosome 
    temp_cna <- cna[cna$chr == chr,]
    temp_snv <- snv[snv$chr == chr,]
    
    # optimized segmentation patameters 
    n_bins <-
      as.integer(chr_info[chr_info$Chromosome == chr,]$Length /
                   fixed_bin_length)
    length_bin <-
      as.integer(chr_info[chr_info$Chromosome == chr,]$Length /
                   n_bins)
    
    end <- length_bin
    start_bin <- NULL
    temp_cna_length <- NULL
    i <- 1
    
    # this function assign each CNA event to segments 
    whichbin <- function(data, end, i) {
      if (any(data$from < end) == TRUE) {
        temp2 <- data[data$from < end,]
        temp2 <-
          cbind(
            temp2,
            start_bin = rep(i, nrow(temp2)),
            bin_length = (temp2$to - temp2$from) / fixed_bin_length
          )
      }
      return(temp2)
    }
    
    # apply the whichbin() function
    for (i in 1:n_bins) {
      if (nrow(temp_cna) == 0) {
        break
      }
      if (any(temp_cna$from < end)) {
        temp_cna_length <-
          rbind(temp_cna_length, whichbin(temp_cna, end, i))
      }
      temp_cna <- temp_cna[!temp_cna$from < end,]
      end <- end + length_bin
    }
    
    # keep CNA events higher that 0,5 Mbp (half of bin size)
    temp_cna_length <-
      temp_cna_length[temp_cna_length$length > fixed_bin_length / 2,]
    amplified_patients <-
      levels(factor(temp_cna_length$ID))
    
    
    ## 2nd loop: load segmented chromosome/gene structure ----
    bin_gene <-
      read.table(
        paste0("/home/ieo5099/Desktop_linux/mutation_compensation/",
               "data/ChromosomeGeneStructure/chr_",
               parse_number(chr),
               "_binSize_",
               fixed_bin_length,
               ".txt"
        )
      )
    
    ## 2nd loop: compute chromosome and arm amplification frequencies ----
    temp_cna_length_backup <- temp_cna_length
    
    ## 2nd loop: compute amplification and deletion frequencies ----
    chr_bins <- data.frame(bins = 1:n_bins)
    
    temp_cna_length <- temp_cna_length_backup
    
    for (i in 1:nrow(temp_cna_length)) {
      if (temp_cna_length$bin_length[i] > 0.5) {
        y <-
          temp_cna_length$bin_length[i] - as.integer(temp_cna_length$bin_length[i])
        if (y >= 0.5) {
          x <-
            data.frame(
              bins = temp_cna_length$start_bin[i]:as.integer(
                temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i]
              ),
              c(rep(
                1, as.integer(temp_cna_length$bin_length[i] + 1)
              ))
            )
        } else{
          x <-
            data.frame(
              bins = temp_cna_length$start_bin[i]:as.integer(
                temp_cna_length$start_bin[i] + temp_cna_length$bin_length[i] - 1
              ),
              c(rep(
                1, as.integer(temp_cna_length$bin_length[i])
              ))
            )
        }
      } else{
        next
      }
      
      chr_bins <- full_join(chr_bins, x, by = "bins")
      
      colnames(chr_bins) <-
        c(colnames(chr_bins)[-length(colnames(chr_bins))], paste0("CNA", i))
    }
    rm(i)
    rm(x)
    rm(y)
    
    chr_bins <- t(as.matrix(chr_bins))
    chr_bins <- as.data.frame(chr_bins)
    chr_bins <- chr_bins[-1,]
    
    # CNA shorter than half of the bin size (0.5 Mbp in this case) were removed from calculation
    chr_bins <-
      cbind(
        total_CN = temp_cna_length[temp_cna_length$bin_length > 0.5,]$total_cn,
        patient = temp_cna_length[temp_cna_length$bin_length > 0.5,]$ID,
        cna_length = temp_cna_length[temp_cna_length$bin_length > 0.5 ,]$length,
        chr_bins
      )
    
    cna_freq_ampl <-
      as.numeric(colSums(chr_bins[chr_bins$total_CN > 2,4:ncol(chr_bins)], na.rm = TRUE) /
                   length(common_patients))
    cna_freq_del <-
      as.numeric(colSums(chr_bins[chr_bins$total_CN < 2,4:ncol(chr_bins)], na.rm = TRUE) /
                   length(common_patients))
    cna_freq_total <-
      as.numeric(colSums(chr_bins[chr_bins$total_CN != 2,4:ncol(chr_bins)], na.rm = TRUE) / length(common_patients))
    
    bin_gene <- cbind(
      bin_gene,
      cna_freq_ampl = cna_freq_ampl,
      cna_freq_del = cna_freq_del,
      cna_freq_total = cna_freq_total
    )
    
    ## 2nd loop: compute the mutation score in copy-neutral regions ----
    # (normalized for the number of patients (1), the coding region (2) and thus log10 (3))
    chr_bins[is.na(chr_bins)] <- 0
    
    chr_bins_pt <- data.frame()
    for (pt in levels(factor(chr_bins$patient))) {
      colSums(chr_bins[chr_bins$patient == pt,][,-c(1:3)])
      chr_bins_pt <-
        rbind(chr_bins_pt, c(pt, as.numeric(colSums(chr_bins[chr_bins$patient == pt,][,-c(1:3)]))))
    }
    colnames(chr_bins_pt) <- c("patients", 1:n_bins)
    
    x <- NULL
    end <- length_bin
    start_bin <- 0
    mutations_raw <- NULL
    mutations_norm <- NULL
    mutations_coding <- NULL
    mutations_noncoding <- NULL
    mutations_coding_norm <- NULL
    mutations_noncoding_norm <- NULL
    n_patients <- NULL
    
    for (i in 1:n_bins) {
      diploid_pt <- length(c(chr_bins_pt[chr_bins_pt[, i + 1] == 0,]$patient))
      mut <-
        temp_snv[temp_snv$ID %in% 
                   c(chr_bins_pt[chr_bins_pt[, i +1] == 0,]$patients),]
      mut <-
        mut[mut$from >= start_bin &
              mut$to < end,]
      
      mut_coding <- nrow(mut[mut$consequence_type %in%
                               names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                            "exon")]) |
                               mut$consequence_type %in%
                               names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                            "missense")]) |
                               mut$consequence_type %in%
                               names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                            "synonymous")]),])
      mut_noncoding <- nrow(mut[mut$consequence_type  %in%
                                  names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                               "intergenic")]) |
                                  mut$consequence_type %in%
                                  names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                               "intron")]) & !(mut$consequence_type %in%
                                                                                                 names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                                                                                              "exon")]) |
                                                                                                 mut$consequence_type %in%
                                                                                                 names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                                                                                              "missense")]) |
                                                                                                 mut$consequence_type %in%
                                                                                                 names(table(mut$consequence_type)[str_detect(names(table(mut$consequence_type)),
                                                                                                                                              "synonymous")])),])

      # mut_coding <- sum(!is.na(mut$transcript_affected))
      # mut_noncoding <- sum(is.na(mut$transcript_affected))
      
      mut_a <- as.data.frame(table(mut$ID))
      
      if (dim(mut_a)[1] != 0) {
        colnames(mut_a) <- c("patient_id", "mutations_raw")
        x <- sum(mut_a$mutations_raw , na.rm = T)
        mutations_raw <- c(mutations_raw, x)
        mutations_norm <- c(mutations_norm, x / diploid_pt)
        mutations_coding <- c(mutations_coding, mut_coding)
        mutations_coding_norm <- c(mutations_coding_norm, mut_coding / diploid_pt)
        mutations_noncoding <- c(mutations_noncoding, mut_noncoding)
        mutations_noncoding_norm <- c(mutations_noncoding_norm, mut_noncoding / diploid_pt)
      } else{
        mutations_raw <- c(mutations_raw, 0)
        mutations_norm <- c(mutations_norm, 0)
        mutations_coding <- c(mutations_coding, 0)
        mutations_coding_norm <- c(mutations_coding_norm, 0)
        mutations_noncoding <- c(mutations_noncoding, 0)
        mutations_noncoding_norm <- c(mutations_noncoding_norm, 0)
        
      }
      n_patients <- c(n_patients, diploid_pt)
      start_bin <- start_bin + length_bin
      end <- end + length_bin
    }
    
    bin_gene_mut <-
      cbind(bin_gene, 
            mutations_raw, 
            mutations_norm, 
            mutations_coding,
            mutations_coding_norm,
            mutations_noncoding,
            mutations_noncoding_norm,
            n_patients)
    bin_gene_mut <- as.data.frame(bin_gene_mut)
    
    bin_gene_mut$length_noncoding <- (bin_gene_mut$bin_end-bin_gene_mut$bin_start)-bin_gene_mut$length_coding
    bin_gene_mut$length_noncoding <- ifelse(bin_gene_mut$length_noncoding < 0, 0, bin_gene_mut$length_noncoding)
    
    ## end of the 2nd loop: write chromosome table ----
    write.table(
      bin_gene_mut,
      file = paste0(
        results_table_path,
        tumor_type,"_",
        chr,
        "_",
        as.integer(fixed_bin_length),
        "BIN_table.txt"
      )
    )
  })
  # }
  print(paste0("\n >>>> ", tumor_type, " DONE! \n"))
}
)

rm(list=ls())
gc(full=T)