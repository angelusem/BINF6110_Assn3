# BINF6110 Assignment 3: draft 

## Introduction
This repository contains files for Assignment 3 in BINF6110. The project analyzes metagenomic sequencing data to compare vegan and omnivore samples.

## Methods 
The overall workflow
1. selection of 3 vegan and 3 omnivore samples
2. quality control with FastQC
3. taxonomic classification with Kraken2
4. abundance estimation with Bracken
5. downstream analysis and visualization in R

Six shotgun metagenome samples from SRA study SRP126540, were selected for this analysis: three omnivore and three vegan subjects. To reduce non-diet confounding, samples were chosen from the same collection site (Turin, Italy) with similar reported sequencing characteristics

Repository structure
- `meta/` : sample metadata and accession information
- `scripts/` : analysis scripts
- `results/figures/` : generated figures
- `results/tables/` : generated tables
- `results/qc/` : FastQC outputs
- `results/kraken2/` : Kraken2 classification results
- `results/bracken/` : Bracken abundance estimates


## Results

Figure 1: Genus-level taxonomic composition
Relative abundance of the most abundant genera across six gut microbiome samples, grouped by diet. Bars show sample-level community composition at the genus level, with low-abundance taxa combined into “Other.” Profiles indicate substantial inter-individual variability, with some diet-associated tendencies but no complete visual separation between groups.

Figure 2: Alpha diversity by diet
Alpha-diversity comparisons between omnivore and vegan samples using Observed richness, Shannon diversity, Simpson diversity, and Berger–Parker dominance. Points represent individual samples and boxplots summarize group distributions. Diversity estimates overlapped substantially between diets, indicating no clear and consistent difference in within-sample diversity.

Figure 3: Bray–Curtis beta-diversity ordination
Principal coordinates analysis of Bray–Curtis dissimilarities at the genus level. Samples show partial separation by diet, suggesting a moderate shift in community composition associated with dietary group, although within-group variability remains evident.

Figure 4: Differential abundance at genus level
ANCOM-BC2 log fold changes for genera comparing vegan versus omnivore samples. Points indicate estimated effect sizes and horizontal bars show confidence intervals. Only g__Segatella was robustly significant, with lower abundance in vegans and higher abundance in omnivores.



## Discussion
Diet appeared to influence gut community composition more than within-sample diversity in this analysis. Alpha-diversity metrics overlapped substantially between vegans and omnivores, with no clear separation in Observed richness, Shannon diversity, or Simpson diversity, while Bray–Curtis PERMANOVA suggested a moderate but non-significant diet effect, indicating that any diet-associated signal was modest relative to inter-individual variation. This cautious interpretation aligns with previous human literature: a systematic review found no consistent microbiome signature separating vegan or vegetarian diets from omnivorous diets across studies (), and a comparative study in urban U.S. vegans and omnivores reported that metabolomic differences were clearer than gut-microbiota differences ().

The clearest signal in the present study was differential abundance of g__Segatella, which was lower in vegans and higher in omnivores by ANCOM-BC2 (lfc_dietvegan = -5.55, q = 0.027). This pattern was supported by the underlying Kraken/Bracken outputs and by the sample-level abundance pattern, where all three omnivore samples had much higher Segatella relative abundance than all three vegan samples. Leave-one-out reanalysis preserved the negative direction of effect in every run, although robust significance was lost after removing any single sample, indicating that the direction of association was stable but statistical certainty was limited by sample size. This finding is biologically plausible given that the original source cohort study () reported diet-associated structuring of Prevotella copri strain repertoires, a bacterium present in the Segatella genus, with omnivore-associated strains enriched for functions related to branched-chain amino acid biosynthesis and vegan-associated strains enriched for genes involved in complex carbohydrate degradation. However, the original study also found that overall P. copri abundance was not significantly associated with diet across the full cohort, emphasizing that the most meaningful differences may occur at the strain and functional level rather than at the genus level alone. Because the present analysis used Kraken2/Bracken genus-level classification on only six samples, the observed Segatella enrichment in omnivores should be interpreted conservatively. No trimming was performed before Kraken2/Bracken classification. This may have added some uncertainty to downstream taxonomic and abundance results, since poor-quality reads can reduce analysis quality and trimming is one way to remove low-quality or contaminating sequence. This limitation may reflect unresolved QC irregularities rather than obvious adapter contamination,  since adapter content passed in all files. Future direction could therefore be replicating the analysis on a larger scale, with a bigger database, trimming and on a species or strain level classification.


## References


