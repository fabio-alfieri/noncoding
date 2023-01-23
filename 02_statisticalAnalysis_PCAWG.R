# SEGMENT-LEVEL ANALYSIS
#
# This script produces as output mutation score and amplification frequency for
# each chromosome (for each tumor type) with different bin sizes and conditions

# rm(list=ls())
gc(full=T)

suppressMessages({
  require(ggplot2)
  require(stringr)
  library(readr)
  require(mgcv)
  require(ggpubr)
  library(optparse)
  library(dplyr)
  library(utils)
})

# set optparse parameters 
option_list = list(
  make_option(
    c("-t", "--tables"),
    type = "character",
    default = "y",
    help = "Options are: ([y]/n)
                It requires some minutes (up to 1 h). 
                It produces bin level correlation tables for different segment
                lengths and conditions.",
    metavar = ""
  ),
  make_option(
    c("-s", "--statistics"),
    type = "character",
    default = "y",
    help = "Options are: ([y]/n)
                It produces statistics and plots.",
    metavar = ""
  )
)


opt_parser = OptionParser(option_list = option_list)

opt = parse_args(opt_parser)


if (opt$tables == "y" & opt$statistics == "y") {
  cat("\n\n >> You chose default options: \n\t(1) --tables 'y';\n\t(2) --statistics 'y' \n\n")
}

if (!any(opt$tables %in% c("y", "n"))) {
  print_help(opt_parser)
  stop("typo in the analysis flag, plase see above for the available options!",
       call. = FALSE)
}

cat(
  "\n\n > This script \n\n\t (1) produce correlations between amplification frequency and mu score for different bin sizes:
    \t\t - from 1 to 50 Mbp
    \t\t - chromosome-arm level
    \t\t - entire chromosome level \n
    \n\t (2) produces correlations and figures \n\n\n"
)

produce_tables <- TRUE
if (opt$tables == "n") {
  produce_tables <- FALSE
}
produce_statistics <- TRUE
if (opt$statistics == "n") {
  produce_statistics <- FALSE
}

# setwd("../")
setwd("~/mountHD/noncoding")

tumor_types <- paste0(c(
  "Breast",
  # "Colorectal",
  "Brain",
  "Kidney",
  "Liver",
  "Pancreas",
  "Lung",
  "Prostate",
  "Ovary",
  "Skin",
  "Blood"
), ".merged")

normLEN <- T

if (F) {
  columns <- read.table("/home/ieo5099/mountHPC/scratch/MS/falfieri/PCAWG/patient_per_tumortype.tsv",
                      skip = 1)
tumor_types <- levels(factor(columns$V1))
tumor_types <- as.data.frame(cbind(do.call(rbind, str_split(tumor_types, "-")), tumor_types))
tumor_types <- tumor_types$tumor_types
tumor_types <- tumor_types[!(tumor_types %in% c("BLCA-US", "BRCA-US", "CESC-US",
                                                "COAD-US", "DLBC-US", "GBM-US",
                                                "HNSC-US", "KICH-US", "KIRC-US",
                                                "KIRP-US", "LAML-KR", "LGG-US", 
                                                "LIHC-US", "LUAD-US", "LUSC-US",
                                                "OV-US", "PRAD-US", "READ-US",
                                                "SARC-US", "SKCM-US", "STAD-US",
                                                "THCA-US", "UCEC-US"
                                                ))]
}

if (produce_tables) {
  # output table directory
  results_table_path <-
    paste0("results/02_produceStatistics/")
  system(paste0("mkdir -p ", results_table_path))
  
  # initialize parameters
  parameters_chr <- data.frame()
  parameters_arm <- data.frame()
  parameters <- data.frame()
  tierfinal <- data.frame()
  
  # segmentation lengths (default is 36Mbp)
  segment_lengths <- c(36)
  
      ## >> 1st loop (optional): varying segment lengths << ----
      for (segment_length in segment_lengths) {
        # 36 was a choose as default (best) cutoff (Fig.1d)
        cat("\n > bin size (Mbp): ", segment_length, "\n")
        
        if (segment_length != 36) {
          conditions <- c(
            "amplifications"
            ,"deletions")
        } else{
          conditions <- c(
            "amplifications"
            ,"deletions"
            ,"coding"
            ,"noncoding"
          )
        }
        
        
        ## >> 2nd loop: mutation/gene type conditions << ----
              source_table_path <-
                paste0("results/01_binLevel/")
          
          ## >> 3rd loop: CANCER TYPE << ----
          suppressMessages({
            for (tumor_type in tumor_types) {
              # 3rd loop: load chromosome tables (output from script 01_mainAnalysis.R) ----
              chr1 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr1_1000000BIN_table.txt"
                  )
                ), chr = 1)
              chr2 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr2_1000000BIN_table.txt"
                  )
                ), chr = 2)
              chr3 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr3_1000000BIN_table.txt"
                  )
                ), chr = 3)
              chr4 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr4_1000000BIN_table.txt"
                  )
                ), chr = 4)
              chr5 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr5_1000000BIN_table.txt"
                  )
                ), chr = 5)
              chr6 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr6_1000000BIN_table.txt"
                  )
                ), chr = 6)
              chr7 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr7_1000000BIN_table.txt"
                  )
                ), chr = 7)
              chr8 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr8_1000000BIN_table.txt"
                  )
                ), chr = 8)
              chr9 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr9_1000000BIN_table.txt"
                  )
                ), chr = 9)
              chr10 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr10_1000000BIN_table.txt"
                  )
                ), chr = 10)
              chr11 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr11_1000000BIN_table.txt"
                  )
                ), chr = 11)
              chr12 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr12_1000000BIN_table.txt"
                  )
                ), chr = 12)
              chr13 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr13_1000000BIN_table.txt"
                  )
                ), chr = 13)
              chr14 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr14_1000000BIN_table.txt"
                  )
                ), chr = 14)
              chr15 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr15_1000000BIN_table.txt"
                  )
                ), chr = 15)
              chr16 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr16_1000000BIN_table.txt"
                  )
                ), chr = 16)
              chr17 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr17_1000000BIN_table.txt"
                  )
                ), chr = 17)
              chr18 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr18_1000000BIN_table.txt"
                  )
                ), chr = 18)
              chr19 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr19_1000000BIN_table.txt"
                  )
                ), chr = 19)
              chr20 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr20_1000000BIN_table.txt"
                  )
                ), chr = 20)
              chr21 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr21_1000000BIN_table.txt"
                  )
                ), chr = 21)
              chr22 <-
                cbind(read.table(
                  file = paste0(
                    source_table_path,
                    tumor_type,
                    "_chr22_1000000BIN_table.txt"
                  )
                ), chr = 22)
              
              
              # 3rd loop: add a new column ----
              # to create segments of desired segment_length and remove bins that do not contains genes
              add.col <- function(df) {
                df <- df[df$gene_count != 0, ]
                df <- cbind(df,
                            resize = rep(
                              1:round(dim(df)[1] / segment_length + 0.5),
                              each = round(dim(df)[1] / round(
                                dim(df)[1] / segment_length + 0.5
                              ) + 0.5)
                            )[1:nrow(df)])
                return(df)
              }
              
              chr1 <-  add.col(chr1)
              chr2 <-  add.col(chr2)
              chr3 <-  add.col(chr3)
              chr4 <-  add.col(chr4)
              chr5 <-  add.col(chr5)
              chr6 <-  add.col(chr6)
              chr7 <-  add.col(chr7)
              chr8 <-  add.col(chr8)
              chr9 <-  add.col(chr9)
              chr10 <- add.col(chr10)
              chr11 <- add.col(chr11)
              chr12 <- add.col(chr12)
              chr13 <- add.col(chr13)
              chr14 <- add.col(chr14)
              chr15 <- add.col(chr15)
              chr16 <- add.col(chr16)
              chr17 <- add.col(chr17)
              chr18 <- add.col(chr18)
              chr19 <- add.col(chr19)
              chr20 <- add.col(chr20)
              chr21 <- add.col(chr21)
              chr22 <- add.col(chr22)
              
              # 3rd loop: create segments of desired segment_length ----
              merge.bins <- function(df) {
                df2 <- data.frame()
                chr <- as.numeric(df$chr[1])
                for (i in levels(factor(df$resize))) {
                  # df <- df[!(df$gene_count <= 1 & df$length_perc <= 0.05), ]
                  
                  start <- df[df$resize == i, ]$bin_start[1]
                  end <-
                    df[df$resize == i, ]$bin_end[length(df[df$resize == i, ]$bin_start)]
                  gene_count <- sum(df[df$resize == i, ]$gene_count)
                  
                  length_coding <- sum(df[df$resize == i, ]$length_coding)
                  length_noncoding <- sum(df[df$resize == i, ]$length_noncoding)
                  
                  cna_freq_ampl <- mean(df[df$resize == i, ]$cna_freq_ampl)
                  cna_freq_del <-  mean(df[df$resize == i, ]$cna_freq_del)
                  
                  mutations_raw <-
                    sum(as.numeric(df[df$resize == i, ]$mutations_raw))
                  mutations_coding <-
                    sum(as.numeric(df[df$resize == i, ]$mutations_coding))
                  mutations_noncoding <-
                    sum(as.numeric(df[df$resize == i, ]$mutations_noncoding))
                  
                  gene_id <-
                    paste0(df[df$resize == i, ]$gene_id, collapse = "")
                  
                  if (dim(df[df$resize == i, ])[1] == 0) {
                    next
                  }
                  if (mutations_raw == 0) {
                    mutations_raw <- 0.0001
                    mutations_coding <- 0.0001
                    mutations_noncoding <- 0.0001
                    # mutations_normPT <- 0.0001
                    # mutations_coding_normPT <- 0.0001
                    # mutations_noncoding_normPT <- 0.0001
                    # mutations_coding_normPT_normLEN <- 0.0001
                    # mutations_noncoding_normPT_normLEN <- 0.0001
                  }
                  mutations_normPT <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_norm),
                      na.rm = T
                    )
                  mutations_coding_normPT <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_coding_norm),
                      na.rm = T
                    )
                  mutations_noncoding_normPT <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_noncoding_norm),
                      na.rm = T
                    )
                  mutations_coding_normPT_normLEN <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_coding_norm) / (as.numeric(df[df$resize == i, ]$length_coding)+1),
                      na.rm = T
                    )
                  mutations_noncoding_normPT_normLEN <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_noncoding_norm) / (as.numeric(df[df$resize == i, ]$length_noncoding)+1),
                      na.rm = T
                    )
                  
                  # mutations are normalized according to Eq. 2
                  
                  df2 <- rbind.data.frame(
                    df2,
                    c(
                      gene_count,
                      start,
                      end,
                      length_coding,
                      length_noncoding,
                      cna_freq_ampl,
                      cna_freq_del,
                      mutations_raw,
                      mutations_coding,
                      mutations_noncoding,
                      mutations_normPT,
                      mutations_coding_normPT,
                      mutations_noncoding_normPT,
                      mutations_coding_normPT_normLEN,
                      mutations_noncoding_normPT_normLEN,
                      gene_id,
                      gene_count,
                      chr,
                      resize = i
                    ),
                    stringsAsFactors = FALSE
                  )
                }
                
                colnames(df2) <- c(
                  "ene_count",
                  "start",
                  "end",
                  "length_coding",
                  "length_noncoding",
                  "cna_freq_ampl",
                  "cna_freq_del",
                  "mutations_raw",
                  "mutations_coding",
                  "mutations_noncoding",
                  "mutations_normPT",
                  "mutations_coding_normPT",
                  "mutations_noncoding_normPT",
                  "mutations_coding_normPT_normLEN",
                  "mutations_noncoding_normPT_normLEN",
                  "gene_id",
                  "gene_count",
                  "chr",
                  "resize"
                )
                
                return(df2)
              }
              
              chr1 <-  merge.bins(chr1)
              chr2 <-  merge.bins(chr2)
              chr3 <-  merge.bins(chr3)
              chr4 <-  merge.bins(chr4)
              chr5 <-  merge.bins(chr5)
              chr6 <-  merge.bins(chr6)
              chr7 <-  merge.bins(chr7)
              chr8 <-  merge.bins(chr8)
              chr9 <-  merge.bins(chr9)
              chr10 <- merge.bins(chr10)
              chr11 <- merge.bins(chr11)
              chr12 <- merge.bins(chr12)
              chr13 <- merge.bins(chr13)
              chr14 <- merge.bins(chr14)
              chr15 <- merge.bins(chr15)
              chr16 <- merge.bins(chr16)
              chr17 <- merge.bins(chr17)
              chr18 <- merge.bins(chr18)
              chr19 <- merge.bins(chr19)
              chr20 <- merge.bins(chr20)
              chr21 <- merge.bins(chr21)
              chr22 <- merge.bins(chr22)
              
              tier <-
                rbind(
                  chr1,
                  chr2,
                  chr3,
                  chr4,
                  chr5,
                  chr6,
                  chr7,
                  chr8,
                  chr9,
                  chr10,
                  chr11,
                  chr12,
                  chr13,
                  chr14,
                  chr15,
                  chr16,
                  chr17,
                  chr18,
                  chr19,
                  chr20,
                  chr21,
                  chr22
                )
              
              tier[, c(1:15, 17:19)] <-
                apply(tier[, c(1:15, 17:19)], 2, as.numeric)
              
              # output table of the analysis
              write_tsv(
                tier,
                file = paste0(
                  results_table_path,
                  tumor_type,
                  "_",
                  segment_length,
                  "Mbp_table.tsv"
                )
              )
              
              
              if(normLEN){
                if(F){
                  corP_amplifications <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                      tier[tier$mutations_normPT != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_amplifications <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                      tier[tier$mutations_normPT != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_deletions <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                 tier[tier$mutations_normPT != 0,]$cna_freq_del, method = "pearson"))
                  corS_deletions <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                 tier[tier$mutations_normPT != 0,]$cna_freq_del, method = "spearman"))
                  
                  corP_coding <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT_normLEN), 
                                              tier[tier$mutations_coding != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_coding <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT_normLEN), 
                                              tier[tier$mutations_coding != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_coding_deletions <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT_normLEN), 
                                                        tier[tier$mutations_coding != 0,]$cna_freq_del, method = "pearson"))
                  corS_coding_deletions <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT_normLEN), 
                                                        tier[tier$mutations_coding != 0,]$cna_freq_del, method = "spearman"))
                  
                  corP_noncoding <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT_normLEN), 
                                                 tier[tier$mutations_noncoding != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_noncoding <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT_normLEN), 
                                                 tier[tier$mutations_noncoding != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_noncoding_deletions <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT_normLEN), 
                                                           tier[tier$mutations_noncoding != 0,]$cna_freq_del, method = "pearson"))
                  corS_noncoding_deletions <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT_normLEN), 
                                                           tier[tier$mutations_noncoding != 0,]$cna_freq_del, method = "spearman"))
                }
                corP_amplifications <- try(cor.test(log10(tier$mutations_normPT), 
                                                    tier$cna_freq_ampl, method = "pearson"))
                corS_amplifications <- try(cor.test(log10(tier$mutations_normPT), 
                                                    tier$cna_freq_ampl, method = "spearman"))
                
                corP_deletions <- try(cor.test(log10(tier$mutations_normPT), 
                                               tier$cna_freq_del, method = "pearson"))
                corS_deletions <- try(cor.test(log10(tier$mutations_normPT), 
                                               tier$cna_freq_del, method = "spearman"))
                
                corP_coding <- try(cor.test(log10(tier$mutations_coding_normPT_normLEN), 
                                            tier$cna_freq_ampl, method = "pearson"))
                corS_coding <- try(cor.test(log10(tier$mutations_coding_normPT_normLEN), 
                                            tier$cna_freq_ampl, method = "spearman"))
                
                corP_coding_deletions <- try(cor.test(log10(tier$mutations_coding_normPT_normLEN), 
                                                      tier$cna_freq_del, method = "pearson"))
                corS_coding_deletions <- try(cor.test(log10(tier$mutations_coding_normPT_normLEN), 
                                                      tier$cna_freq_del, method = "spearman"))
                
                corP_noncoding <- try(cor.test(log10(tier$mutations_noncoding_normPT_normLEN), 
                                               tier$cna_freq_ampl, method = "pearson"))
                corS_noncoding <- try(cor.test(log10(tier$mutations_noncoding_normPT_normLEN), 
                                               tier$cna_freq_ampl, method = "spearman"))
                
                corP_noncoding_deletions <- try(cor.test(log10(tier$mutations_noncoding_normPT_normLEN), 
                                                         tier$cna_freq_del, method = "pearson"))
                corS_noncoding_deletions <- try(cor.test(log10(tier$mutations_noncoding_normPT_normLEN), 
                                                         tier$cna_freq_del, method = "spearman"))
              }else{
                if(F){
                  corP_amplifications <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                      tier[tier$mutations_normPT != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_amplifications <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                      tier[tier$mutations_normPT != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_deletions <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                 tier[tier$mutations_normPT != 0,]$cna_freq_del, method = "pearson"))
                  corS_deletions <- try(cor.test(log10(tier[tier$mutations_normPT != 0,]$mutations_normPT), 
                                                 tier[tier$mutations_normPT != 0,]$cna_freq_del, method = "spearman"))
                  
                  corP_coding <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT), 
                                              tier[tier$mutations_coding != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_coding <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT), 
                                              tier[tier$mutations_coding != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_coding_deletions <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT), 
                                                        tier[tier$mutations_coding != 0,]$cna_freq_del, method = "pearson"))
                  corS_coding_deletions <- try(cor.test(log10(tier[tier$mutations_coding != 0,]$mutations_coding_normPT), 
                                                        tier[tier$mutations_coding != 0,]$cna_freq_del, method = "spearman"))
                  
                  corP_noncoding <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT), 
                                                 tier[tier$mutations_noncoding != 0,]$cna_freq_ampl, method = "pearson"))
                  corS_noncoding <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT), 
                                                 tier[tier$mutations_noncoding != 0,]$cna_freq_ampl, method = "spearman"))
                  
                  corP_noncoding_deletions <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT), 
                                                           tier[tier$mutations_noncoding != 0,]$cna_freq_del, method = "pearson"))
                  corS_noncoding_deletions <- try(cor.test(log10(tier[tier$mutations_noncoding != 0,]$mutations_noncoding_normPT), 
                                                           tier[tier$mutations_noncoding != 0,]$cna_freq_del, method = "spearman"))
                }
                corP_amplifications <- try(cor.test(log10(tier$mutations_normPT), 
                                                    tier$cna_freq_ampl, method = "pearson"))
                corS_amplifications <- try(cor.test(log10(tier$mutations_normPT), 
                                                    tier$cna_freq_ampl, method = "spearman"))
                
                corP_deletions <- try(cor.test(log10(tier$mutations_normPT), 
                                               tier$cna_freq_del, method = "pearson"))
                corS_deletions <- try(cor.test(log10(tier$mutations_normPT), 
                                               tier$cna_freq_del, method = "spearman"))
                
                corP_coding <- try(cor.test(log10(tier$mutations_coding_normPT), 
                                            tier$cna_freq_ampl, method = "pearson"))
                corS_coding <- try(cor.test(log10(tier$mutations_coding_normPT), 
                                            tier$cna_freq_ampl, method = "spearman"))
                
                corP_coding_deletions <- try(cor.test(log10(tier$mutations_coding_normPT), 
                                                      tier$cna_freq_del, method = "pearson"))
                corS_coding_deletions <- try(cor.test(log10(tier$mutations_coding_normPT), 
                                                      tier$cna_freq_del, method = "spearman"))
                
                corP_noncoding <- try(cor.test(log10(tier$mutations_noncoding_normPT), 
                                               tier$cna_freq_ampl, method = "pearson"))
                corS_noncoding <- try(cor.test(log10(tier$mutations_noncoding_normPT), 
                                               tier$cna_freq_ampl, method = "spearman"))
                
                corP_noncoding_deletions <- try(cor.test(log10(tier$mutations_noncoding_normPT), 
                                                         tier$cna_freq_del, method = "pearson"))
                corS_noncoding_deletions <- try(cor.test(log10(tier$mutations_noncoding_normPT), 
                                                         tier$cna_freq_del, method = "spearman"))
              }
              
              tryCatch({
                parameters <- rbind(parameters,  rbind(c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "amplifications",
                  corP = corP_amplifications$estimate,
                  p.corP = corP_amplifications$p.value,
                  corS = corS_amplifications$estimate,
                  p.corS = corS_amplifications$p.value
                ),c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "deletions",
                  corP = corP_deletions$estimate,
                  p.corP = corP_deletions$p.value,
                  corS = corS_deletions$estimate,
                  p.corS = corS_deletions$p.value
                ),c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "coding",
                  corP = corP_coding$estimate,
                  p.corP = corP_coding$p.value,
                  corS = corS_coding$estimate,
                  p.corS = corS_coding$p.value
                ),c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "coding_deletions",
                  corP = corP_coding_deletions$estimate,
                  p.corP = corP_coding_deletions$p.value,
                  corS = corS_coding_deletions$estimate,
                  p.corS = corS_coding_deletions$p.value
                ),c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "noncoding",
                  corP = corP_noncoding$estimate,
                  p.corP = corP_noncoding$p.value,
                  corS = corS_noncoding$estimate,
                  p.corS = corS_noncoding$p.value
                ),c(
                  tumor_type = tumor_type,
                  segment_length = segment_length,
                  condition = "noncoding_deletions",
                  corP = corP_noncoding_deletions$estimate,
                  p.corP = corP_noncoding_deletions$p.value,
                  corS = corS_noncoding_deletions$estimate,
                  p.corS = corS_noncoding_deletions$p.value
                )))
              })
             
              
          }
          })
        
      }
      
      colnames(parameters) <-
        c("tumorType",
          "segment_length",
          "condition",
          "corP",
          "p.corP",
          "corS",
          "p.corS")
      parameters[, 4:7] <- apply(parameters[, 4:7], 2, as.numeric)
      
}
