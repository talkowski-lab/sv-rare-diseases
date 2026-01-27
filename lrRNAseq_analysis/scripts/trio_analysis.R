library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(glue)
library(stringr)

# Table of transcript_id\tgene_id\tgene_name for gencode 48
transcript_table <- ""
id_table <- read_delim(transcript_table) %>%
  dplyr::select(-transcript_id) %>%
  distinct() %>%
  mutate(gene_id_simple = sapply(strsplit(gene_id, ".", fixed=TRUE), function(s) s[1]))




variants <- read_delim('../data/lrRNA_Trios/denovo_snvs_lrRNAseq_samples.txt') 
families <- variants %>% 
  pull(SAMPLE) %>%
  unique() %>%
  gsub("_3$", "", .) 

all_exp <- lapply(families, function(fam) {
  lapply(1:3, function(i) {

    indiv <- glue("{fam}_{i}")
    isoclass <- read_delim(glue('../results/isoform_match/sqanti3_qc/{indiv}/{indiv}_classification.txt'))
    exp <- read_delim(glue('../data/lrRNA_Trios/transcriptomes/{indiv}_R1.LRAA.quant.expr')) %>%
      left_join(
        isoclass %>% select(isoform, structural_category, associated_gene, associated_transcript, predicted_NMD),
        by = c('transcript_id' = 'isoform')
      ) %>%
      left_join(
        id_table, by=c('associated_gene' = 'gene_id')
      ) %>%
      mutate(indiv = indiv, family = fam, relation = case_when(
        i == 1 ~ "Parent1",
        i == 2 ~ "Parent2",
        i == 3 ~ "Child"
      ))
    exp
  }) %>% bind_rows()
}) %>% bind_rows()

all_junc <- lapply(families, function(fam) {
  lapply(1:3, function(i) {
    indiv <- glue("{fam}_{i}")
    juncclass <- read_delim(glue('../results/isoform_match/sqanti3_qc/{indiv}/{indiv}_junctions.txt')) %>%
      mutate(indiv = indiv, family=fam)
  }) %>% bind_rows()
}) %>% bind_rows()


all_ml_filter <- lapply(families, function(fam) {
  lapply(1:3, function(i) {
    indiv <- glue("{fam}_{i}")
    read_delim(glue("../results/isoform_match/sqanti3_qc/{indiv}/filter_ML/{indiv}_MLresult_classification.txt")) %>%
      select(isoform, filter_result) %>%
      mutate(indiv = indiv, family=fam) 
  }) %>% bind_rows()
}) %>% bind_rows()

all_basic_filter <- lapply(families, function(fam) {
  lapply(1:3, function(i) {
    indiv <- glue("{fam}_{i}")
    read_delim(glue("../results/isoform_match/sqanti3_qc/{indiv}/filter_rules/{indiv}_RulesFilter_result_classification.txt")) %>%
      select(isoform, filter_result) %>%
      mutate(indiv = indiv, family=fam) 
  }) %>% bind_rows()
}) %>% bind_rows()


all_exp <- all_exp %>% 
  left_join(all_ml_filter %>% filter(filter_result == "Isoform") %>% select(transcript_id=isoform, indiv) %>% mutate(filter_ML=TRUE)) %>%
  left_join(all_basic_filter %>% filter(filter_result == "Isoform") %>% select(transcript_id=isoform, indiv) %>% mutate(filter_rules=TRUE)) %>%
  replace_na(list(filter_ML=FALSE, filter_rules=FALSE))

all_exp <- all_exp %>%
  group_by(indiv, associated_gene) %>%
  mutate(gene_TPM = sum(TPM)) %>%
  ungroup()

median_gene_exp <- all_exp %>%
  select(gene_TPM, indiv, associated_gene) %>%
  distinct() %>%
  group_by(associated_gene) %>%
  summarize(median_gene_TPM = median(gene_TPM))


fusions <- all_exp %>% 
  filter(structural_category == "fusion") %>%
  separate_longer_delim(associated_gene, delim="_") %>%
  select(-gene_name) %>%
  left_join(id_table, by=c('associated_gene' = 'gene_id')) 




fusions %>% 
  filter(filter_ML, filter_rules, relation == "Child")

fusions %>%
  select(transcript_id, indiv) %>%
  distinct() %>%
  group_by(indiv) %>%
  tally() %>%
  mutate(dummy = 0) %>%
  ggplot(aes(x=n)) +
    geom_blank(aes(x=dummy)) +
    geom_histogram(binwidth=1, color="black") +
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
    theme_clean(base_size=14) 
ggsave("../results/isoform_match/figures/fusion_transcript_counts.pdf", width=3, height=3)

fusions %>%
  filter(filter_rules, filter_ML) %>%
  select(transcript_id, indiv) %>%
  distinct() %>%
  group_by(indiv) %>%
  tally() %>%
  mutate(dummy = 0) %>%
  ggplot(aes(x=n)) +
    geom_blank(aes(x=dummy)) +
    geom_histogram(binwidth=1, color="black") +
    scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
    theme_clean(base_size=14) 
ggsave("../results/isoform_match/figures/fusion_filtered_transcript_counts.pdf", width=3, height=3)


# Make some QC figures on number of transcripts found
high_exp_genes <- median_gene_exp %>%
  filter(median_gene_TPM > 1) %>%
  pull(associated_gene) %>%
  unique()



categs <- c("full-splice_match",  "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic_intron", "genic")
all_categs <- c(categs, "antisense", "fusion", "intergenic")
categ_colors <- c("#6baed6", "#fc8d59", "#78c679", "#ee6a50", "#a427ec", "#969696")

all_categ_colors <- c(categ_colors, "#fb52a1", "#b1a192", "#212121")

mean_frac <- lapply(categs, function(categ) {
  all_exp %>%
    group_by(indiv, associated_gene) %>%
    summarize(has_categ = any(structural_category == categ)) %>%
    group_by(indiv) %>%
    summarize(
      n_genes = n(),
      n_genes_with_categ = sum(has_categ),
      frac_genes_with_categ = n_genes_with_categ / n_genes
    ) %>%
    mutate(structural_category = categ)
}) %>% bind_rows()

mean_frac_hq <- lapply(categs, function(categ) {
  all_exp %>%
    filter(associated_gene %in% high_exp_genes) %>%
    group_by(indiv, associated_gene) %>%
    summarize(has_categ = any(structural_category == categ & unique_gene_read_fraction > 0.1 & filter_ML & filter_rules)) %>%
    group_by(indiv) %>%
    summarize(
      n_genes = n(),
      n_genes_with_categ = sum(has_categ),
      frac_genes_with_categ = n_genes_with_categ / n_genes
    ) %>%
    mutate(structural_category = categ)
}) %>% bind_rows()

mean_frac %>%
  mutate(dummy=0) %>%
  mutate(structural_category = factor(structural_category, levels=categs)) %>%
  ggplot(aes(x=n_genes_with_categ)) +
  geom_histogram(aes(fill=structural_category)) +
  geom_blank(aes(x=dummy)) +
  scale_fill_manual(values=categ_colors) +
  facet_wrap(structural_category ~ ., ncol=1) +
  theme_clean(base_size=14) + 
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )
ggsave("../results/isoform_match/figures/frac_structural_categs.pdf", width=3, height=6)

mean_frac_hq %>%
  mutate(dummy=0) %>%
  mutate(structural_category = factor(structural_category, levels=categs)) %>%
  ggplot(aes(x=n_genes_with_categ)) +
  geom_histogram(aes(fill=structural_category)) +
  geom_blank(aes(x=dummy)) +
  scale_fill_manual(values=categ_colors) +
  facet_wrap(structural_category ~ ., ncol=1, scales="free_x") +
  theme_clean(base_size=14) + 
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )
ggsave("../results/isoform_match/figures/frac_structural_categs_hq.pdf", width=3, height=6)

all_exp %>%
  filter(!is.na(structural_category)) %>%
  group_by(indiv, structural_category) %>%
  tally() %>%
  mutate(dummy=0) %>%
  mutate(structural_category = factor(structural_category, levels=all_categs)) %>%
  ggplot(aes(x=n)) +
  geom_histogram(aes(fill=structural_category)) +
  geom_blank(aes(x=dummy)) +
  scale_fill_manual(values=all_categ_colors) +
  labs(x="# of Transcripts", y="# of Samples") +
  facet_wrap(structural_category ~ ., scales="free", labeller=labeller(structural_category=setNames(str_to_title(gsub("[-_]", " ", all_categs)), all_categs))) +
  theme_clean(base_size=14) + 
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )
ggsave("../results/isoform_match/figures/n_transcript_types_hist.pdf", width=8, height=6)

all_exp %>%
  filter(filter_ML, filter_rules) %>%
  filter(!is.na(structural_category)) %>%
  group_by(indiv, structural_category) %>%
  tally() %>%
  mutate(dummy=0) %>%
  mutate(structural_category = factor(structural_category, levels=all_categs)) %>%
  ggplot(aes(x=n)) +
  geom_histogram(aes(fill=structural_category)) +
  geom_blank(aes(x=dummy)) +
  scale_fill_manual(values=all_categ_colors) +
  labs(x="# of Transcripts", y="# of Samples") +
  facet_wrap(structural_category ~ ., scales="free", labeller=labeller(structural_category=setNames(str_to_title(gsub("[-_]", " ", all_categs)), all_categs))) +
  theme_clean(base_size=14) + 
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle=30, hjust=1),
    panel.spacing = unit(0, "lines"),
    legend.position = "none"
  )
ggsave("../results/isoform_match/figures/n_transcript_types_hist_filtered.pdf", width=8, height=6)

# Get all transcripts within genes with de novo variants
dd <- all_exp %>% 
  filter(relation == "Child") %>%
  mutate(gene_id_simple = sapply(strsplit(associated_gene, ".", fixed=TRUE), function(s) s[1])) %>%
  left_join(variants %>%
    select(Gene, ID, SAMPLE) %>% filter(!is.na(Gene)),
    by=c('gene_id_simple' = 'Gene', 'indiv' = 'SAMPLE'),
    relationship="many-to-many"
  ) %>%
  filter(!is.na(ID))


novel_dd <- dd %>% # Novel transcripts in genes with de novo variants
  filter(associated_transcript == "novel")


dd_fusions <- fusions %>% # Fusion transcripts in genes with de novo variants
  filter(relation == "Child") %>%
  mutate(gene_id_simple = sapply(strsplit(associated_gene, ".", fixed=TRUE), function(s) s[1])) %>%
  left_join(variants %>%
    select(Gene, ID, SAMPLE) %>% filter(!is.na(Gene)),
    by=c('gene_id_simple' = 'Gene', 'indiv' = 'SAMPLE'),
    relationship="many-to-many"
  ) %>%
  filter(!is.na(ID))

dd_fusions %>% 
  select(indiv, gene_name, transcript_id, structural_category, TPM, gene_TPM, ID) %>%
  write_delim("../results/denovo_fusion_transcripts.tsv", delim="\t")


# Filter out transcripts that are in parents
novel_dd <- lapply(split(novel_dd, novel_dd$family), function(df) {
  fam <- unique(df$family)
  parents <- all_exp %>%
    filter(family == fam, relation != "Child") 
  df %>%
    filter(!(exons %in% parents$exons))
}) %>% bind_rows()

length(variants$ID)
variants %>% filter(!is.na(Feature)) %>% nrow()
sum(variants$ID %in% novel_dd$ID)
get_indivs_and_nvar <- function(df) {
  c(
    df %>% nrow(),
    df %>% pull(SAMPLE) %>% unique() %>% length()
  )
}
get_indivs_and_nvar(variants %>% filter(ID %in% novel_dd$ID))
get_indivs_and_nvar(variants %>% filter(ID %in% (novel_dd %>% filter(TPM > 1) %>% pull(ID))))
get_indivs_and_nvar(variants %>% filter(ID %in% (novel_dd %>% filter(TPM > 1, isoform_fraction > 0.1) %>% pull(ID))))

get_indivs_and_nvar(variants %>% filter(ID %in% (novel_dd %>% filter(filter_ML, filter_rules) %>% pull(ID))))
get_indivs_and_nvar(variants %>% filter(ID %in% (novel_dd %>% filter(filter_ML, filter_rules, TPM > 1) %>% pull(ID))))
get_indivs_and_nvar(variants %>% filter(ID %in% (novel_dd %>% filter(filter_ML, filter_rules, TPM > 1, isoform_fraction > 0.1) %>% pull(ID))))



novel_dd %>%
  filter(filter_ML, filter_rules, isoform_fraction > 0.1) %>%
  select(indiv, gene_name, transcript_id, structural_category, TPM, gene_TPM, ID) %>%
  write_delim("../results/denovo_novel_transcripts.tsv", delim="\t")



### Now examine novel transcripts near inherited recessive variants with at least one hit
rare_var_vep <- read_delim("../data/lrRNA_Trios/rare_variants_vep_highmod.txt") %>%
  rename(ID = `#Uploaded_variation`) %>%
  select(ID, gene_id_simple=Gene, Consequence, Extra) %>%
  distinct() %>%
  mutate(indiv = sapply(strsplit(ID, ":", fixed=TRUE), function(s) s[5])) %>%
  mutate(family = sapply(strsplit(indiv, "_", fixed=TRUE), function(s) paste(c(s[1], s[2]), collapse="_"))) %>%
  left_join(id_table %>% select(gene_id_simple, gene_name))

rare_var_vep_high <- rare_var_vep %>%
  filter(grepl("HIGH", Extra))

rd <- all_exp %>%
  filter(relation == "Child") %>%
  mutate(gene_id_simple = sapply(strsplit(associated_gene, ".", fixed=TRUE), function(s) s[1])) %>%
  left_join(rare_var_vep_high %>%
    select(gene_id_simple, ID, indiv, Consequence),
    by=c('gene_id_simple', 'indiv'),
    relationship="many-to-many"
  ) %>%
  filter(!is.na(ID))


novel_rd <- rd %>%
  filter(associated_transcript == "novel")

rd_fusions <- fusions %>%
  filter(relation == "Child") %>%
  mutate(gene_id_simple = sapply(strsplit(associated_gene, ".", fixed=TRUE), function(s) s[1])) %>%
  left_join(rare_var_vep_high %>%
    select(gene_id_simple, ID, indiv, Consequence),
    by=c('gene_id_simple', 'indiv'),
    relationship="many-to-many"
  ) %>%
  filter(!is.na(ID))

rd_fusions %>% 
  select(indiv, gene_name, transcript_id, structural_category, TPM, gene_TPM, ID) %>%
  write_delim("../results/rarevar_fusion_transcripts.tsv", delim="\t")
rd_vars <- rare_var_vep_high  %>%
  select(ID, indiv) %>%
  distinct()

length(rd_vars$ID)
sum(rd_vars$ID %in% novel_rd$ID)
get_indivs_and_nvar <- function(df) {
  c(
    df %>% nrow(),
    df %>% pull(indiv) %>% unique() %>% length()
  )
}
get_indivs_and_nvar(rd_vars %>% filter(ID %in% novel_rd$ID))
get_indivs_and_nvar(rd_vars %>% filter(ID %in% (novel_rd %>% filter(TPM > 1) %>% pull(ID))))
get_indivs_and_nvar(rd_vars %>% filter(ID %in% (novel_rd %>% filter(TPM > 1, isoform_fraction > 0.1) %>% pull(ID))))

get_indivs_and_nvar(rd_vars %>% filter(ID %in% (novel_rd %>% filter(filter_ML, filter_rules) %>% pull(ID))))
get_indivs_and_nvar(rd_vars %>% filter(ID %in% (novel_rd %>% filter(filter_ML, filter_rules, TPM > 1) %>% pull(ID))))
get_indivs_and_nvar(rd_vars %>% filter(ID %in% (novel_rd %>% filter(filter_ML, filter_rules, TPM > 1, isoform_fraction > 0.1) %>% pull(ID))))

novel_rd %>%
  filter(filter_ML, filter_rules, isoform_fraction > 0.1) %>%
  select(indiv, gene_name, transcript_id, structural_category, TPM, gene_TPM, ID) %>%
  distinct() %>%
  write_delim("../results/rarevar_novel_transcripts.tsv", delim="\t")

# Get recessive variants close to splice sites
interval_dist <- function(a,b, x,y) {
  d1 <- x - b
  d2 <- y - a
  if (sign(d1) != sign(d2)) {
    return(0)  # intervals overlap
  } else {
    return(min(abs(d1), abs(d2)))
  }
}
recessive_splice_sites <- lapply(families, function(family) {
  proband <- glue("{family}_3")
  read_delim(glue("../results/isoform_match/sqanti3_qc/{proband}/{proband}_junc_recessive_variant_intersect.bed"), col_names=c("_chrom", "_start", "_end", "junc_id", "ID",  "wgs_start", "wgs_end" ), delim="\t") %>%
    mutate(distance_to_site = mapply(interval_dist, (`_start` + `_end`)/2-1, (`_start` + `_end`)/2+1, wgs_start, wgs_end)) %>%
    select(!(starts_with("_"))) %>%
    mutate(
      isoform = sapply(strsplit(junc_id, ":", fixed=TRUE), function(s) paste(s[1:5], collapse=":")),
      junction_number = sapply(strsplit(junc_id, ":", fixed=TRUE), function(s) paste0("junction_", s[6])),
    ) %>%
    left_join(all_junc %>% filter(indiv == proband), by=c("isoform", "junction_number")) %>%
    left_join(all_exp %>% filter(indiv == proband) %>% select(isoform=transcript_id, associated_gene, associated_transcript), by=c("isoform"))
}) %>% bind_rows()


second_hit_sites <- recessive_splice_sites %>%
  filter(!grepl("novel", associated_gene)) %>%
  mutate(is_fusion = grepl("_", associated_gene)) %>%
  separate_longer_delim(associated_gene, delim="_") %>%
  left_join(id_table, by=c('associated_gene' = 'gene_id')) %>%
  filter(!(ID %in% rd_vars$ID)) %>% # I'm looking for second hits, so exclude the variants in the recessive list
  inner_join(rare_var_vep_high %>% select(lof_ID=ID, gene_name, indiv, Consequence) %>% distinct())


second_hit_sites %>%
  select(indiv, gene_name, junc_id, is_fusion, query_ID=ID, distance_to_site, pLoF_ID=lof_ID) %>%
  write_delim("../results/recessive_splice_sites_second_hit_candidates.tsv", delim="\t")
### Investigate specific list
outrider_results <- read_delim("../results/long_read_gene_level__ods__q2_all_results.tsv")

examine_df <- read_delim("../data/examine_list.txt")
# Create IGV sessions

if (file.exists("../results/igv_sessions/template.xml")) {
  # If there is a template to build igv sessions to investigate variants
  template <- readLines("../results/igv_sessions/template.xml")
  for (i in 1:nrow(examine_df)) {
    family_id <- gsub("RGP_", "", examine_df$family[i], fixed=TRUE)
    if (file.exists(glue("../results/igv_sessions/{family_id}_{examine_df$gene_name[i]}.xml"))) next
    gene_name <- examine_df$gene_name[i]
    sapply(template, glue) %>% writeLines(glue("../results/igv_sessions/{family_id}_{gene_name}.xml"))
  }
}

results_df <- lapply(1:nrow(examine_df), function(i) {
  gene  <- examine_df$gene_id[i]
  gene_name <- examine_df$gene_name[i]
  family <- examine_df$family[i]
  individual <- glue("{family}_3")


  filtered_exp <- all_exp %>% filter(
    gene_name == !!gene_name,
    family == !!family
  )
  child_exp <- filtered_exp %>% filter(relation == "Child")
  p1_exp <- filtered_exp %>% filter(relation == "Parent1")
  p2_exp <- filtered_exp %>% filter(relation == "Parent2")

  # Get OUTRIDER z-score and p-value
  zsc <- outrider_results %>%
    filter(sampleID == individual, symbol == gene_name) %>%
    pull(zScore)
  pv <- outrider_results %>%
    filter(sampleID == individual, symbol == gene_name) %>%
    pull(pValue)
  zsc <- ifelse(length(zsc) == 0, NA, zsc)
  pv <- ifelse(length(pv) == 0, NA, pv)

  # Get TPM
  child_tpm <- child_exp %>% pull(TPM) %>% sum()
  p1_tpm <- p1_exp %>% pull(TPM) %>% sum()
  p2_tpm <- p2_exp %>% pull(TPM) %>% sum()

  # Get fusion transcripts
  child_fusion <- all_exp %>% filter(indiv == individual, structural_category == "fusion") %>% filter(grepl(gene, associated_gene)) %>% pull(transcript_id)
  child_fusion <- ifelse(length(child_fusion) == 0, NA, child_fusion)

  # Get fraction non full splice match
  child_frac <- child_exp %>% filter(structural_category != "full-splice_match") %>% pull(isoform_fraction) %>% sum()
  p1_frac <- p1_exp %>% filter(structural_category != "full-splice_match") %>% pull(isoform_fraction) %>% sum()
  p2_frac <- p2_exp %>% filter(structural_category != "full-splice_match") %>% pull(isoform_fraction) %>% sum()

  # Same but filtered
  child_frac_filt <- (child_exp %>% filter(structural_category != "full-splice_match", filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum()) /
  (child_exp %>% filter(filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum())
  p1_frac_filt <- (p1_exp %>% filter(structural_category != "full-splice_match", filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum()) /
  (p1_exp %>% filter(filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum())
  p2_frac_filt <- (p2_exp %>% filter(structural_category != "full-splice_match", filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum()) /
  (p2_exp %>% filter(filter_ML, filter_rules) %>% pull(isoform_fraction) %>% sum())

  # Get any high expressed non full splice match isoforms
  child_high_exp_iso <- child_exp %>%
    filter(structural_category != "full-splice_match") %>%
    filter(TPM > 0.5, isoform_fraction > 0.10) %>%
    pull(transcript_id)
  child_high_exp_frac <- child_exp %>%
    filter(structural_category != "full-splice_match") %>%
    filter(TPM > 0.5, isoform_fraction > 0.10) %>%
    pull(isoform_fraction) 
  child_high_exp_frac <- round(child_high_exp_frac * 100, 1)



  p1_has_high_exp_iso <- sapply(child_high_exp_iso, function(iso) {
    exons <- child_exp %>% filter(transcript_id == iso) %>% pull(exons)
    return(exons %in% p1_exp$exons)
  })

  p2_has_high_exp_iso <- sapply(child_high_exp_iso, function(iso) {
    exons <- child_exp %>% filter(transcript_id == iso) %>% pull(exons)
    return(exons %in% p2_exp$exons)
  })

  child_high_exp_iso <- if_else(length(child_high_exp_iso) == 0, NA_character_,
    paste(glue("{child_high_exp_iso} ({child_high_exp_frac}%)"), collapse=", ")
  )

  p1_has_high_exp_iso <- if_else(length(p1_has_high_exp_iso) == 0, NA_character_, paste(p1_has_high_exp_iso, collapse=", "))
  p2_has_high_exp_iso <- if_else(length(p2_has_high_exp_iso) == 0, NA_character_, paste(p2_has_high_exp_iso, collapse=", "))

  # Get any recessive variant overlap
  child_junc_rec <- recessive_splice_sites %>% 
    filter(indiv == individual, grepl(gene, associated_gene)) %>%
    pull(junc_id)
  child_junc_rec <- ifelse(length(child_junc_rec) == 0, NA, paste(child_junc_rec, collapse=", "))

  setNames(c(gene, gene_name, individual, family, 
    zsc, pv,
    child_tpm, p1_tpm, p2_tpm,
    child_fusion, 
    child_frac, p1_frac, p2_frac,
    child_frac_filt, p1_frac_filt, p2_frac_filt,
    child_high_exp_iso, p1_has_high_exp_iso, p2_has_high_exp_iso,
    child_junc_rec),
  c("gene_id", "gene_name", "individual", "family",
    "outrider_zscore", "outrider_pvalue",
    "child_TPM", "parent1_TPM", "parent2_TPM",
    "child_fusion_transcripts",
    "child_frac_non_full_splice_match", "parent1_frac_non_full_splice_match", "parent2_frac_non_full_splice_match",
    "child_frac_non_full_splice_match_filt", "parent1_frac_non_full_splice_match_filt", "parent2_frac_non_full_splice_match_filt",
    "child_high_exp_isoforms", "p1_has_high_exp_isoforms", "p2_has_high_exp_isoforms",
    "recessive_junctions_overlapping_with_gene")
  )
}) %>% bind_rows()
results_df %>%
  write_delim("../results/gene_check.tsv", delim="\t")


