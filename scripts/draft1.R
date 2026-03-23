# ASSIGNMENT 3 : BINF6110

# packages 
library(biomformat)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tibble)
library(vegan)
library(ANCOMBC)

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

# only 9 taxa across the 6 samples, but proceeding for now
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

wilcox.test(Shannon ~ diet, data = alpha_df)
wilcox.test(Observed ~ diet, data = alpha_df)