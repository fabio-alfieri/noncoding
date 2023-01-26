corParameters <- parameters
corParameters$tumorType <- as.factor(corParameters$tumorType)
corParameters$condition <- as.factor(corParameters$condition)
                                  
corParameters$tumorType <- factor(corParameters$tumorType,
                                  levels = c(corParameters[corParameters$condition == "coding", ][order(corParameters[corParameters$condition == "coding", ]$corS, decreasing = T), ]$tumorType))


# figure synonymous vs. non_synonymous (revisions) ----
synonym <-
  corParameters[
    # corParameters$condition == "amplifications" |
                  # corParameters$condition == "deletions" |
                  # corParameters$condition == "coding_deletions" |
                  # corParameters$condition == "noncoding_deletions" 
                  corParameters$condition == "coding" |
                  # corParameters$condition == "coding_wointron" |
                  # corParameters$condition == "noncoding_intron" |
                  corParameters$condition == "noncoding_intergenic"
                  # corParameters$condition == "deletions"
                  ,]

p1 <-
  ggplot(synonym, aes(
    x = as.factor(tumorType),
    y = corS,
    fill = as.factor(condition)
  )) +
  geom_bar(#color = "black",
    stat = "identity",
    position = position_dodge()) +
  theme_classic() + scale_fill_brewer(palette = "Blues") +
  theme(legend.position = "top")
p2 <-
  ggplot(synonym, aes(
    x = as.factor(tumorType),
    y = log10(p.corS),
    fill = condition
  )) +
  geom_bar(#color = "black",
    stat = "identity",
    position = position_dodge()) +
  theme_classic() + scale_fill_brewer(palette = "Reds") +
  geom_hline(yintercept = log10(0.1)) +
  theme(legend.position = "bottom")

# ggarrange(p1, p2, nrow = 2, ncol = 1)

  p3 <- ggplot(synonym, aes(
    x = corS,
    y = condition,
    fill = as.factor(condition))) +
    geom_boxplot() +
    theme_classic() +
    scale_fill_brewer(palette = "Blues") +
    geom_jitter() + coord_flip() +
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


print(ggarrange(ggarrange(p1, p2, nrow = 2, ncol = 1), 
                ggarrange(NULL, p3, NULL, ncol = 1, nrow = 3, heights = c(0.3, 1, 0.3)), 
                widths = c(1, 0.2)))


wilcox.test(synonym[synonym$condition == "coding",]$corS,
            synonym[synonym$condition == "noncoding",]$corS, paired = T)

