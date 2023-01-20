
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

tumor_types <- c(
  "PACA-CA"
)

segment_cutoffs <- c(20)

stringent_mutations <- F
if(stringent_mutations){
  segment_cutoffs_muts <- c( #2.5,
    # 5,
    10,
    15,
    20
    # 30
  )
}else{
  segment_cutoffs_muts <- 20
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
  
  for(segment_cutoff_muts in segment_cutoffs_muts){
    cat("\n > MUTATION segment_mean cutoff: ", segment_cutoff_muts)
    
    for(segment_cutoff in segment_cutoffs){
      if(stringent_mutations){
        cat("\n > COPY-NUMBER segment_mean cutoff: ", segment_cutoff, "\n")
      }
      cat("\n >> (1) analyis: produce correlations at different bin sizes \n")
      
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
          )
        }
        
        if (segment_length == 1) {
          # produce chromosome/arm level correlations
          chr_arm <- "yes" # "" or "yes"
        } else{
          chr_arm <- ""
        }
        
        ## >> 2nd loop: mutation/gene type conditions << ----
        for (condition in conditions) {
              source_table_path <-
                paste0("results/01_binLevel/")
          
          cat(" > correlation condition: ", condition, "\n")
          
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
                  df <- df[!(df$gene_count <= 1 & df$length_perc <= 0.05), ]
                  start <- df[df$resize == i, ]$bin_start[1]
                  end <-
                    df[df$resize == i, ]$bin_end[length(df[df$resize == i, ]$bin_start)]
                  gene_count <- sum(df[df$resize == i, ]$gene_count)
                  length_coding <-
                    sum(df[df$resize == i, ]$length_coding)
                  length_perc <- mean(df[df$resize == i, ]$length_perc)
                  cna_freq_ampl <-
                    mean(df[df$resize == i, ]$cna_freq_ampl)
                  cna_freq_del <- mean(df[df$resize == i, ]$cna_freq_del)
                  mutations_raw <-
                    sum(as.numeric(df[df$resize == i, ]$mutations_raw))
                  gene_id <-
                    paste0(df[df$resize == i, ]$gene_id, collapse = "")
                  if (mutations_raw == 0) {
                    mutations_raw <- 0
                  }
                  if (dim(df[df$resize == i, ])[1] == 0) {
                    next
                  }
                  mutations_norm <-
                    mean(
                      as.numeric(df[df$resize == i, ]$mutations_norm) / as.numeric(df[df$resize == i, ]$length_coding),
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
                      length_perc,
                      cna_freq_ampl,
                      cna_freq_del,
                      mutations_raw,
                      mutations_norm,
                      gene_id,
                      chr,
                      resize = i
                    ),
                    stringsAsFactors = FALSE
                  )
                }
                
                colnames(df2) <- c(
                  "gene_count",
                  "start",
                  "end",
                  "length_coding",
                  "length_perc",
                  "cna_freq_ampl",
                  "cna_freq_del",
                  "mutations_raw",
                  "mutations_norm",
                  "gene_id",
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
              
              tier[, c(1:9, 11:12)] <-
                apply(tier[, c(1:9, 11:12)], 2, as.numeric)
              
              # output table of the analysis
              write_tsv(
                tier,
                file = paste0(
                  results_table_path,
                  condition,
                  "_",
                  tumor_type,
                  "_",
                  segment_length,
                  "Mbp_table.tsv"
                )
              )
              
              # produce plots
              if (condition != "deletions") {
                x <- log10(tier[tier$mutations_raw != 0, ]$mutations_raw)
                y <- tier[tier$mutations_norm != 0, ]$cna_freq_ampl
                
                corP <- try(cor.test(x, y, method = "pearson"))
                corS <- try(cor.test(x, y, method = "spearman"))
                
                parameters <- rbind(
                  parameters,
                  c(
                    tumor_type,
                    segment_length,
                    paste0(condition),
                    corP$estimate,
                    p.corP = corP$p.value,
                    corS$estimate,
                    p.corS = corS$p.value
                  ),
                  stringsAsFactors = FALSE
                )
                (
                  p1 <- ggplot(tier[tier$mutations_norm != 0, ], aes(x = x, y = y)) +
                    geom_point() +
                    theme_classic() +
                    geom_smooth(method = "lm", alpha = 0.15) +
                    # geom_smooth(method = "gam", formula = y ~ s(x), alpha = 0.05, se = FALSE, size = 1, color = "black") +
                    ggtitle(
                      paste(tumor_type, "-", segment_length, "-", condition),
                      subtitle = paste(
                        "Pearson's R =",
                        signif(corP$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corP$p.value, digits = 3),
                        "\nSpearman's rho =",
                        signif(corS$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corS$p.value, digits = 3)
                      )
                    ) +
                    xlab("Log of diploid mutations") +
                    ylab("Amplification frequency")
                )
                
                system(
                  paste0(
                    "mkdir -p results/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/"
                  )
                )
                pdf(
                  file = paste0(
                    "results/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/",
                    tumor_type,
                    "_",
                    condition,
                    ".pdf"
                  )
                )
                print(p1)
                dev.off()
              }
              
              if (condition == "deletions") {
                x <- log10(tier[tier$mutations_raw != 0, ]$mutations_raw)
                y <- tier[tier$mutations_norm != 0, ]$cna_freq_del
                
                corP <- try(cor.test(x, y, method = "pearson"))
                corS <- try(cor.test(x, y, method = "spearman"))
                
                parameters <- rbind(
                  parameters,
                  c(
                    tumor_type,
                    segment_length,
                    paste0(condition),
                    corP$estimate,
                    p.corP = corP$p.value,
                    corS$estimate,
                    p.corS = corS$p.value
                  ),
                  stringsAsFactors = FALSE
                )
                (
                  p1 <- ggplot(tier[tier$mutations_norm != 0, ], aes(x = x, y = y)) +
                    geom_point() +
                    theme_classic() +
                    geom_smooth(method = "lm", alpha = 0.15) +
                    # geom_smooth(method = "gam", formula = y ~ s(x), alpha = 0.05, se = FALSE, size = 1, color = "black") +
                    ggtitle(
                      paste(tumor_type, "-", segment_length, "-", condition),
                      subtitle = paste(
                        "Pearson's R =",
                        signif(corP$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corP$p.value, digits = 3),
                        "\nSpearman's rho =",
                        signif(corS$estimate, digits = 3),
                        "\n",
                        "p-value =",
                        signif(corS$p.value, digits = 3)
                      )
                    ) +
                    xlab("Log of diploid mutations") +
                    ylab("Deletion frequency")
                )
                
                system(
                  paste0(
                    "mkdir -p results/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/"
                  )
                )
                pdf(
                  file = paste0(
                    "results/02_tumorCorrelations/",
                    segment_length,
                    "Mbp_",
                    condition,
                    "/",
                    tumor_type,
                    "_",
                    condition,
                    ".pdf"
                  )
                )
                print(p1)
                dev.off()
              }
            }
          })
        }
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
  }
}
