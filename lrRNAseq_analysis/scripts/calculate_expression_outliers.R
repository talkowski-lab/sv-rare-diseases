library(OUTRIDER)
library(BiocParallel)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggthemes)
library(glue)


outrider_results <- read_delim("../results/long_read_gene_level__ods__q2_all_results.tsv")

gene_df <- outrider_results %>%
  select(geneID, symbol) %>%
  distinct() %>%
  group_by(symbol) %>%
  filter(n() == 1) %>%
  ungroup()

gene_dict <- setNames(gene_df$geneID, gene_df$symbol)

zsc_mat <- outrider_results %>% 
  select(sampleID, geneID, zScore) %>%
  group_by(sampleID, geneID, zScore) %>%
  filter(n() == 1) %>%
  pivot_wider(names_from="sampleID", values_from="zScore") %>%
  arrange(geneID) %>%
  column_to_rownames("geneID") %>%
  as.matrix()

pval_mat <- outrider_results %>% 
  select(sampleID, geneID, pValue) %>%
  group_by(sampleID, geneID, pValue) %>%
  filter(n() == 1) %>%
  pivot_wider(names_from="sampleID", values_from="pValue") %>%
  arrange(geneID) %>%
  column_to_rownames("geneID") %>%
  as.matrix()



maz_id <- gene_df %>% filter(symbol == "MAZ") %>% pull(geneID)

data.frame(z=zsc_mat[maz_id, ], indiv=colnames(zsc_mat)) %>%
  ggplot(aes(x = reorder(indiv, z), y = z)) +
  geom_hline(yintercept=0, color="grey") +
  geom_point(aes(color = (indiv == "RGP_2159_3"))) +
  scale_color_manual(values=c("black", "red")) +
  labs(x="Individual", y="Z-score") +
  theme_clean() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="none"
  )
ggsave("../results/isoform_match/figures/MAZ_expression.pdf", width=3, height=3)


dmr_gene_names <- readLines("../data/gene_lists/all_DMR_genes.txt")
dmr_gene_ids <- gene_dict[dmr_gene_names]

zsc_indiv <- zsc_mat[, "RGP_2159_3"]

zsc_dmr <- zsc_indiv[dmr_gene_ids]
mean(zsc_dmr, na.rm=TRUE)
median(zsc_dmr)
sd(zsc_dmr) / sqrt(length(zsc_dmr))

set.seed(620)
N <- 10000
length <- sum(!is.na(zsc_dmr))
samples <- replicate(N, {
  sample_zsc <- sample(zsc_indiv, length, replace=TRUE)
  c(mean(sample_zsc, trim=0.1), sd(sample_zsc), median(sample_zsc), mean(sample_zsc^2), median(sample_zsc^2))
})

sim_means <- samples[1,]
sim_sd <- samples[2,]
sim_medians <- samples[3,]
sim_means2 <- samples[4,]
sim_medians2 <- samples[5,]


mean(zsc_dmr, na.rm=TRUE, trim=0.1)
sd(zsc_dmr, na.rm=TRUE)

p_mean <- mean(sim_means > mean(zsc_dmr, na.rm=TRUE), na.rm=TRUE)
p_sd <- mean(sim_sd > sd(zsc_dmr, na.rm=TRUE), na.rm=TRUE)
p_medians <- mean(sim_medians > median(zsc_dmr, na.rm=TRUE), na.rm=TRUE)
p_mean2 <- mean(sim_means2 > mean(zsc_dmr ^ 2, na.rm=TRUE), na.rm=TRUE)
p_median2 <- mean(sim_medians2 > median(zsc_dmr ^ 2, na.rm=TRUE), na.rm=TRUE)

ggplot(data.frame(x=zsc_dmr), aes(x=x)) +
  geom_histogram() +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="Z-score") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_DMR_zscore_histogram.pdf", width=3, height=1.8)

ggplot(data.frame(x=zsc_indiv), aes(x=x)) +
  geom_histogram() +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="Z-score") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_all_zscore_histogram.pdf", width=3, height=1.8)


ggplot(data.frame(x=sim_means), aes(x=x)) +
  geom_histogram() +
  geom_vline(xintercept=mean(zsc_dmr, na.rm=TRUE), color="red", linewidth=2) +
  geom_text(label = glue("p={p_mean}"), x = - 0.05, y = 1000, color="red") + 
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="Mean Z-score") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_meanDMR_zscore_histogram.pdf", width=3, height=1.8)

ggplot(data.frame(x=sim_medians), aes(x=x)) +
  geom_histogram(binwidth=0.01) +
  geom_vline(xintercept=median(zsc_dmr, na.rm=TRUE), color="red", linewidth=2) +
  geom_text(label = glue("p={p_medians}"), x = - 0.1, y = 1000, color="red") + 
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="Median Z-score") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_medianDMR_zscore_histogram.pdf", width=3, height=1.8)

data.frame(sd=sim_sd) %>%
  ggplot(aes(x=sd)) +
  geom_histogram() +
  geom_vline(xintercept=sd(zsc_dmr, na.rm=TRUE), color="red", linewidth=2) +
  geom_text(label = glue("p={p_sd}"), x = 1.075, y = 1250, color="red") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="SD Z-Score") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_sdDMR_zscore_histogram.pdf", width=3, height=1.8)

data.frame(mean2=sim_means2) %>%
  ggplot(aes(x=mean2)) +
  geom_histogram() +
  geom_vline(xintercept=mean(zsc_dmr^2, na.rm=TRUE), color="red", linewidth=2) +
  geom_text(label = glue("p={p_mean2}"), x = 0.99, y = 900, color="red") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="Mean(Z-Score^2)") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_meanDMR_zscoresquared_histogram.pdf", width=3, height=1.8)

data.frame(median2=sim_medians2) %>%
  ggplot(aes(x=median2)) +
  geom_histogram(binwidth=0.01) +
  geom_vline(xintercept=median(zsc_dmr^2, na.rm=TRUE), color="red", linewidth=2) +
  geom_text(label = glue("p={p_median2}"), x = 0.62, y = 1250, color="red") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  labs(x="median(Z-Score^2)") +
  theme_clean()
ggsave("../results/isoform_match/figures/RGP_2159_3_medianDMR_zscoresquared_histogram.pdf", width=3, height=1.8)
