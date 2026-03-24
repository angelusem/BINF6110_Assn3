# ASSIGNMENT 3 : BINF6110

# packages 
library(biomformat)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tibble)
library(vegan)
library(ANCOMBC)
library(patchwork)

# importing BIOM

biom_data <- read_biom("results/biom/genus_FtoG_json.biom")
physeq <- import_biom(biom_data)

# importing metadata
meta <- read.delim("meta/metadata.tsv", row.names = 1, check.names = FALSE)

# keeping only shared samples, and putting metadata in the same order as phyloseq samples
shared <- intersect(sample_names(physeq), rownames(meta))
physeq <- prune_samples(shared, physeq)
meta <- meta[shared, , drop = FALSE]
sample_data(physeq) <- sample_data(meta)

# exploring this object
physeq
sample_names(physeq)
sample_variables(physeq)
rank_names(physeq)
head(as(sample_data(physeq), "data.frame"))

# 103 taxa across the 6 samples
#
tax_df <- as.data.frame(tax_table(physeq))
head(tax_df, 10)
table(is.na(tax_df$Rank6))

#
# quick checks
ntaxa(physeq)
nsamples(physeq)
sample_sums(physeq)
sort(taxa_sums(physeq), decreasing = TRUE)

# inspecting the taxonomy table
as.data.frame(tax_table(physeq))

# creating output folders
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# seting diet as an ordered factor
sample_data(physeq)$diet <- factor(sample_data(physeq)$diet,
                                   levels = c("omnivore", "vegan"))
# TAXONOMIC ABUNDANCE: Figure 1

physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_gen <- tax_glom(physeq_rel, taxrank = "Rank6", NArm = FALSE)

df_gen <- psmelt(physeq_gen)
df_gen$Rank6 <- as.character(df_gen$Rank6)
df_gen$Rank6[is.na(df_gen$Rank6) | df_gen$Rank6 == ""] <- "Unclassified"

top_taxa <- df_gen %>%
  group_by(Rank6) %>%
  summarise(total_abund = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = 8) %>%
  pull(Rank6)
top_taxa

df_gen$Genus_plot <- ifelse(df_gen$Rank6 %in% top_taxa, df_gen$Rank6, "Other")

p_tax <- ggplot(df_gen, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ diet, scales = "free_x") +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_tax

ggsave("results/figures/taxonomic_abundance_rank6.png", p_tax,
       width = 10, height = 6, dpi = 300)

write.csv(df_gen, "results/tables/taxonomic_abundance_long_rank6.csv", row.names = FALSE)

# For the compuation of Alpha Diversity
alpha_df <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson")) %>%
  rownames_to_column("Sample") %>%
  left_join(as(sample_data(physeq), "data.frame") %>% rownames_to_column("Sample"),
            by = "Sample")
# Adding Berger-Parker dominance 
otu_mat <- as(otu_table(physeq), "matrix")

# Ensure rows = samples, columns = taxa
if (taxa_are_rows(physeq)) {
  otu_mat <- t(otu_mat)
}

berger_parker <- apply(otu_mat, 1, function(x) {
  x <- as.numeric(x)
  if (sum(x) == 0) return(NA_real_)
  max(x) / sum(x)
})

alpha_df$Berger_Parker <- berger_parker[alpha_df$Sample]
alpha_df

write.csv(alpha_df, "results/tables/alpha_diversity_values.csv", row.names = FALSE)

p_alpha_shannon <- ggplot(alpha_df, aes(x = diet, y = Shannon, color = diet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3) +
  theme_bw() +
  labs(x = "Diet", y = "Shannon diversity")

p_alpha_shannon
ggsave("results/figures/alpha_diversity_shannon.png", p_alpha_shannon,
       width = 6, height = 5, dpi = 300)

p_alpha_observed <- ggplot(alpha_df, aes(x = diet, y = Observed, color = diet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3) +
  theme_bw() +
  labs(x = "Diet", y = "Observed richness")

p_alpha_observed
ggsave("results/figures/alpha_diversity_observed.png", p_alpha_observed,
       width = 6, height = 5, dpi = 300)

p_alpha_bp <- ggplot(alpha_df, aes(x = diet, y = Berger_Parker, color = diet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3) +
  theme_bw() +
  labs(x = "Diet", y = "Berger-Parker dominance")

p_alpha_bp
ggsave("results/figures/alpha_diversity_berger_parker.png", p_alpha_bp,
       width = 6, height = 5, dpi = 300)

p_alpha_simpson <- ggplot(alpha_df, aes(x = diet, y = Simpson, color = diet)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 3) +
  theme_bw() +
  labs(x = "Diet", y = "Simpson diversity")

p_alpha_simpson
ggsave("results/figures/alpha_diversity_simpson.png", p_alpha_simpson,
       width = 6, height = 5, dpi = 300)

# saving a summary table
alpha_summary <- alpha_df |>
  dplyr::group_by(diet) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_observed = mean(Observed),
    sd_observed = sd(Observed),
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    mean_simpson = mean(Simpson),
    sd_simpson = sd(Simpson),
    mean_berger_parker = mean(Berger_Parker),
    sd_berger_parker = sd(Berger_Parker)
  )

alpha_summary
write.csv(alpha_summary, "results/tables/alpha_diversity_summary_by_diet.csv", row.names = FALSE)
# running group tests for alpha
wilcox.test(Shannon ~ diet, data = alpha_df)
wilcox.test(Observed ~ diet, data = alpha_df)
wilcox.test(Simpson ~ diet, data = alpha_df)
wilcox.test(Berger_Parker ~ diet, data = alpha_df)

# combining the figures for alpha diversity
p_alpha_observed <- p_alpha_observed + labs(title = "Observed richness")
p_alpha_shannon  <- p_alpha_shannon  + labs(title = "Shannon diversity")
p_alpha_simpson  <- p_alpha_simpson  + labs(title = "Simpson diversity")
p_alpha_bp       <- p_alpha_bp       + labs(title = "Berger-Parker dominance")

alpha_panel <- ((p_alpha_observed | p_alpha_shannon) /
                  (p_alpha_simpson | p_alpha_bp)) +
  plot_annotation(tag_levels = "A")

alpha_panel

ggsave("results/figures/alpha_diversity_panel.png", alpha_panel,
       width = 10, height = 8, dpi = 300)



# Beta diversity and PERMANOVA
ord_bray <- ordinate(physeq, method = "PCoA", distance = "bray")

p_beta <- plot_ordination(physeq, ord_bray, color = "diet") +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "Bray-Curtis PCoA")

p_beta
ggsave("results/figures/beta_diversity_bray_pcoa.png", p_beta,
       width = 7, height = 5, dpi = 300)

dist_bray <- phyloseq::distance(physeq, method = "bray")
meta_df <- as(sample_data(physeq), "data.frame")

set.seed(6110)
perm <- adonis2(dist_bray ~ diet, data = meta_df)

perm
capture.output(perm, file = "results/tables/permanova_bray.txt")
# Jaccard as well
ord_jacc <- ordinate(physeq, method = "PCoA", distance = "jaccard")

p_beta_jacc <- plot_ordination(physeq, ord_jacc, color = "diet") +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "Jaccard PCoA")

p_beta_jacc
ggsave("results/figures/beta_diversity_jaccard_pcoa.png", p_beta_jacc,
       width = 7, height = 5, dpi = 300)

# Differential Abundance Analysis!

# the reference level must be omnivore
sample_data(physeq)$diet <- relevel(factor(sample_data(physeq)$diet), ref = "omnivore")

ancombc.out <- ancombc2(
  data = physeq,
  tax_level = "Rank6",
  fix_formula = "diet",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "diet",
  struc_zero = TRUE,
  neg_lb = TRUE
)

names(ancombc.out)
names(ancombc.out$res)
head(ancombc.out$res)
head(ancombc.out$zero_ind)

# saving the results table
da_res <- as.data.frame(ancombc.out$res)

write.csv(da_res, "results/tables/ancombc2_full_results.csv", row.names = FALSE)

head(da_res)

# to extract significant taxa
sig_da <- da_res |>
  dplyr::filter(diff_robust_dietvegan == TRUE) |>
  dplyr::arrange(q_dietvegan)

sig_da


write.csv(sig_da, "results/tables/ancombc2_significant_taxa.csv", row.names = FALSE)

da_res |>
  dplyr::filter(q_dietvegan < 0.05) |>
  dplyr::arrange(q_dietvegan)

# for the differential abundance plot
p_da <- ggplot(da_res, aes(x = lfc_dietvegan, y = reorder(taxon, lfc_dietvegan))) +
  geom_point(aes(color = diff_robust_dietvegan), size = 3) +
  geom_errorbar(aes(xmin = lfc_dietvegan - se_dietvegan,
                    xmax = lfc_dietvegan + se_dietvegan),
                width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Log fold change (vegan vs omnivore)",
       y = "Taxon",
       color = "Robustly significant")

p_da

ggsave("results/figures/differential_abundance_ancombc2.png", p_da,
       width = 8, height = 7, dpi = 300)
# --- end of analysis ----