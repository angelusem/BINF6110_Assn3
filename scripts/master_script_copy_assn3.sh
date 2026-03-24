# ASSIGNMENT 3: BINF6110

# scripts submitted on Narval:
# 1. downloading reads
cat > $SCRATCH/assn3_binf6110/scripts/fasterq.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=a3_fasterq
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=$SCRATCH/assn3_binf6110/logs/fasterq_%j.out
#SBATCH --error=$SCRATCH/assn3_binf6110/logs/fasterq_%j.err

module load sra-toolkit/3.0.9

cd $SCRATCH/assn3_binf6110
mkdir -p fastq tmp logs

while read ACC; do
    echo "Converting $ACC"
    fasterq-dump raw/$ACC \
        --split-files \
        --threads 8 \
        --temp tmp \
        --outdir fastq
done < meta/SraAccList.txt
EOF

# missing accension from first job, so one-job script for each missing accension

cat > $SCRATCH/assn3_binf6110/scripts/fasterq_one.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=a3_fq_one
#SBATCH --account=def-cottenie
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/%u/assn3_binf6110/logs/fasterq_%j.out
#SBATCH --error=/scratch/%u/assn3_binf6110/logs/fasterq_%j.err

module load sra-toolkit/3.0.9

cd $SCRATCH/assn3_binf6110
mkdir -p fastq tmp logs

rm -f fastq/${ACC}_1.fastq fastq/${ACC}_2.fastq fastq/${ACC}.fastq

echo "Converting ${ACC}"
fasterq-dump raw/${ACC} \
    --split-files \
    --threads 4 \
    --temp tmp \
    --outdir fastq
EOF

# submitting using the list of missed accensions
while read ACC; do
    sbatch --export=ALL,ACC=$ACC $SCRATCH/assn3_binf6110/scripts/fasterq_one.slurm
done < $SCRATCH/assn3_binf6110/meta/missing_acc.txt

# to monitor them
squeue -u $USER

# to verify if I have all 12 fastq files once completed:
find $SCRATCH/assn3_binf6110/fastq -type f | wc -l

for acc in SRR8146971 SRR8146969 SRR8146976 SRR8146974 SRR8146973 SRR8146963; do
    ls -lh $SCRATCH/assn3_binf6110/fastq/${acc}_*.fastq
done
#
# 2. fastqc check 
cat > $SCRATCH/assn3_binf6110/scripts/qc_fastqc_only.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=a3_fastqc
#SBATCH --account=def-cottenie
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/scratch/%u/assn3_binf6110/logs/fastqc_%j.out
#SBATCH --error=/scratch/%u/assn3_binf6110/logs/fastqc_%j.err

set -euo pipefail

module purge
module load StdEnv/2023 fastqc/0.12.1

BASE=$SCRATCH/assn3_binf6110
mkdir -p $BASE/results/qc/fastqc
mkdir -p $BASE/logs

fastqc \
  --threads $SLURM_CPUS_PER_TASK \
  --outdir $BASE/results/qc/fastqc \
  $BASE/fastq/*.fastq
EOF
# submitting job
sbatch $SCRATCH/assn3_binf6110/scripts/qc_fastqc_only.slurm

# checking whether fully finished
ls $SCRATCH/assn3_binf6110/results/qc/fastqc/*_fastqc.html | wc -l
ls $SCRATCH/assn3_binf6110/results/qc/fastqc/*_fastqc.zip | wc -l

# for a summary table
mkdir -p $SCRATCH/assn3_binf6110/results/qc/fastqc_summaries

for z in $SCRATCH/assn3_binf6110/results/qc/fastqc/*_fastqc.zip; do
    bn=$(basename "${z%.zip}")
    unzip -p "$z" "${bn}/summary.txt" > \
      $SCRATCH/assn3_binf6110/results/qc/fastqc_summaries/${bn}_summary.txt
done

cat $SCRATCH/assn3_binf6110/results/qc/fastqc_summaries/*_summary.txt > \
  $SCRATCH/assn3_binf6110/results/qc/fastqc_summaries/all_fastqc_summary.txt

#
# 3. Kracken2 + Bracken analysis

# small submission to diagnose any database issues if any
cd $SCRATCH/assn3_binf6110

cat > scripts/kraken_bracken_test.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=k2test
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/k2test_%j.out
#SBATCH --error=logs/k2test_%j.err

set -euo pipefail

cd $SCRATCH/assn3_binf6110

module purge
module load StdEnv/2023 gcc/12.3 kraken2/2.1.6 bracken/3.0

DB=$SCRATCH/db/k2_standard_08_GB
READ_LEN=150
ACC=SRR8146963

mkdir -p results/kraken2
mkdir -p results/bracken/species
mkdir -p results/bracken/genus

kraken2 \
  --db "$DB" \
  --paired \
  --confidence 0.15 \
  --memory-mapping \
  --threads "$SLURM_CPUS_PER_TASK" \
  --output "results/kraken2/${ACC}.kraken" \
  --report "results/kraken2/${ACC}.kreport" \
  "fastq/${ACC}_1.fastq" "fastq/${ACC}_2.fastq"

bracken \
  -d "$DB" \
  -i "results/kraken2/${ACC}.kreport" \
  -o "results/bracken/species/${ACC}.species.bracken" \
  -r "$READ_LEN" \
  -l S

bracken \
  -d "$DB" \
  -i "results/kraken2/${ACC}.kreport" \
  -o "results/bracken/genus/${ACC}.genus.bracken" \
  -r "$READ_LEN" \
  -l G
EOF

sbatch scripts/kraken_bracken_test.slurm

# test timed out in 4 hours so trying again
cd $SCRATCH/assn3_binf6110

cat > scripts/kraken_bracken_test.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=k2test
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/k2test_%j.out
#SBATCH --error=logs/k2test_%j.err

set -euo pipefail

cd $SCRATCH/assn3_binf6110

module --force purge
module load StdEnv/2023 gcc/12.3 kraken2/2.1.6 bracken/3.0

DB=$SCRATCH/db/k2_standard_08_GB
READ_LEN=150
ACC=SRR8146963

mkdir -p results/kraken2
mkdir -p results/bracken/species
mkdir -p results/bracken/genus

rm -f "results/kraken2/${ACC}.kraken"
rm -f "results/kraken2/${ACC}.kreport"
rm -f "results/bracken/species/${ACC}.species.bracken"
rm -f "results/bracken/genus/${ACC}.genus.bracken"

kraken2 \
  --db "$DB" \
  --paired \
  --confidence 0.15 \
  --memory-mapping \
  --threads "$SLURM_CPUS_PER_TASK" \
  --output "results/kraken2/${ACC}.kraken" \
  --report "results/kraken2/${ACC}.kreport" \
  "fastq/${ACC}_1.fastq" "fastq/${ACC}_2.fastq"

bracken \
  -d "$DB" \
  -i "results/kraken2/${ACC}.kreport" \
  -o "results/bracken/species/${ACC}.species.bracken" \
  -r "$READ_LEN" \
  -l S

bracken \
  -d "$DB" \
  -i "results/kraken2/${ACC}.kreport" \
  -o "results/bracken/genus/${ACC}.genus.bracken" \
  -r "$READ_LEN" \
  -l G
EOF

sbatch scripts/kraken_bracken_test.slurm

# running all 6 files as an array job: 
cat > $SCRATCH/assn3_binf6110/scripts/kraken2_bracken_array.slurm <<'EOF'
#!/bin/bash
#SBATCH --job-name=a3_k2b
#SBATCH --account=def-cottenie
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=1-6
#SBATCH --output=/scratch/%u/assn3_binf6110/logs/k2b_%A_%a.out
#SBATCH --error=/scratch/%u/assn3_binf6110/logs/k2b_%A_%a.err

set -euo pipefail

module purge
module load kraken2 bracken   # replace with exact module names if needed

BASE=$SCRATCH/assn3_binf6110
DB=$SCRATCH/db/k2_standard_08_GB
READ_LEN=150

mkdir -p $BASE/results/kraken2
mkdir -p $BASE/results/bracken/species
mkdir -p $BASE/results/bracken/genus

mapfile -t ACCS < $BASE/meta/SraAccList.txt
ACC=${ACCS[$SLURM_ARRAY_TASK_ID-1]}

kraken2 \
  --db $DB \
  --paired \
  --confidence 0.15 \
  --memory-mapping \
  --threads $SLURM_CPUS_PER_TASK \
  --output $BASE/results/kraken2/${ACC}.kraken \
  --report $BASE/results/kraken2/${ACC}.kreport \
  $BASE/fastq/${ACC}_1.fastq $BASE/fastq/${ACC}_2.fastq

bracken \
  -d $DB \
  -i $BASE/results/kraken2/${ACC}.kreport \
  -o $BASE/results/bracken/species/${ACC}.species.bracken \
  -w $BASE/results/bracken/species/${ACC}.species.kreport \
  -r $READ_LEN \
  -l S

bracken \
  -d $DB \
  -i $BASE/results/kraken2/${ACC}.kreport \
  -o $BASE/results/bracken/genus/${ACC}.genus.bracken \
  -w $BASE/results/bracken/genus/${ACC}.genus.kreport \
  -r $READ_LEN \
  -l G

echo "Finished ${ACC}"
EOF

sbatch $SCRATCH/assn3_binf6110/scripts/kraken2_bracken_array.slurm

# verify outputs: I need to see 6, 6 and 6 to be satisfied and move on to next step
ls $SCRATCH/assn3_binf6110/results/kraken2/*.kreport | wc -l
ls $SCRATCH/assn3_binf6110/results/bracken/species/*.bracken | wc -l
ls $SCRATCH/assn3_binf6110/results/bracken/genus/*.bracken | wc -l

# now that kracken2 and bracken have run succesfully, the next step is working on abundance tables. using the BIOM route

# converting my original csv metadata file to tsv since biom requires tsv:
tr ',' '\t' < meta/metadata.csv > meta/metadata.tsv

# moving all my genus .kreport files to one location 
mkdir -p results/biom_input/genus

ln -s ../../bracken/genus/SRR8146969.genus.kreport results/biom_input/genus/SRR8146969.kreport
ln -s ../../bracken/genus/SRR8146971.genus.kreport results/biom_input/genus/SRR8146971.kreport
ln -s ../../bracken/genus/SRR8146973.genus.kreport results/biom_input/genus/SRR8146973.kreport
ln -s ../../bracken/genus/SRR8146974.genus.kreport results/biom_input/genus/SRR8146974.kreport
ln -s ../../bracken/genus/SRR8146976.genus.kreport results/biom_input/genus/SRR8146976.kreport
ln -s ../../kraken2/SRR8146963_bracken_genuses.kreport results/biom_input/genus/SRR8146963.kreport

ls -l results/biom_input/genus

# for the kracken-biom install, setting up a python Virtual Environment in home directory NOT scratch
module load python
virtualenv --no-download ~/venvs/krakenbiom
source ~/venvs/krakenbiom/bin/activate

# making sure it is in home 
 echo "$VIRTUAL_ENV"
readlink -f "$VIRTUAL_ENV"

# installing kracken-biom
pip install --upgrade pip
pip install kraken-biom
kraken-biom -h

# building the BIOM file, using the json format despite the h5py package being available for the readibility of json format and compatibility with R version 
mkdir -p results/biom

kraken-biom results/biom_input/genus/*.kreport \
  --fmt json \
  -o results/biom/genus_FtoG_json.biom

'''
# default hdf5 format 
mkdir -p results/biom
kraken-biom results/biom_input/genus/*.kreport \
  -o results/biom/genus_table.biom
'''

# confirming
ls -lh results/biom/genus_table.biom

# importing into R
setwd("C:/Users/angel/binf6110/assignments/BINF6110_Assn3")
getwd()
list.files()

# confirming the Segatella result 
# 1) Confirming the direction directly from the abundance table not just in ancombc2
# Relative abundance already exists in df_gen from your script
seg_rel <- df_gen %>%
  dplyr::filter(Rank6 == "g__Segatella") %>%
  dplyr::select(Sample, diet, Abundance) %>%
  dplyr::arrange(diet, desc(Abundance))

seg_rel

# checking also the raw counts at genus level instead of just proportions
physeq_gen_counts <- tax_glom(physeq, taxrank = "Rank6", NArm = FALSE)

tax_tab <- as.data.frame(tax_table(physeq_gen_counts))
otu_mat <- as(otu_table(physeq_gen_counts), "matrix")
if (taxa_are_rows(physeq_gen_counts)) {
  otu_mat <- t(otu_mat)
}

seg_taxa <- rownames(tax_tab)[tax_tab$Rank6 == "g__Segatella"]

seg_counts <- data.frame(
  Sample = rownames(otu_mat),
  Segatella_count = rowSums(otu_mat[, seg_taxa, drop = FALSE]),
  Total_reads = rowSums(otu_mat),
  diet = as(sample_data(physeq_gen_counts), "data.frame")[rownames(otu_mat), "diet"]
)

seg_counts$Segatella_rel <- seg_counts$Segatella_count / seg_counts$Total_reads
seg_counts[order(seg_counts$diet, -seg_counts$Segatella_rel), ]

# for numeric confirmation that omnivore samples mostly higher than vegan samples, and
#not just one tiny count difference,
#not a case where everything is driven by one sample with the others all near zero

# 2) Making sure the taxon is coming from the correct final object not the 9 taxa one
# is the signal present in the final physeq object?
ntaxa(physeq)
nsamples(physeq)
sample_names(physeq)
sample_data(physeq)
table(as.data.frame(tax_table(physeq))$Rank6 == "g__Segatella", useNA = "ifany")
# are sample labels aligned as expected
data.frame(
  Sample = sample_names(physeq),
  diet = sample_data(physeq)$diet,
  site = sample_data(physeq)$site,
  subject_id = sample_data(physeq)$subject_id
)
# sanity checking the direction on ancombc2
wilcox.test(Segatella_rel ~ diet, data = seg_counts)
tapply(seg_counts$Segatella_rel, seg_counts$diet, summary)

# to check whether one sample is driving everything

# on a plot
ggplot(seg_counts, aes(x = diet, y = Segatella_rel, label = Sample, color = diet)) +
  geom_point(size = 3) +
  geom_text(nudge_x = 0.1, size = 3, show.legend = FALSE) +
  theme_bw()

# more thorough l-o-o method
loo_results <- lapply(sample_names(physeq), function(s) {
  phy_sub <- prune_samples(sample_names(physeq) != s, physeq)
  sample_data(phy_sub)$diet <- relevel(factor(sample_data(phy_sub)$diet), ref = "omnivore")

  out <- ancombc2(
    data = phy_sub,
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

  res <- as.data.frame(out$res)
  seg <- res[res$taxon == "g__Segatella", c("taxon", "lfc_dietvegan", "q_dietvegan", "diff_robust_dietvegan")]
  seg$left_out <- s
  seg
})

doo <- do.call(rbind, loo_results)
doo

# checking on source files (kraken2 and bracken)
grep -H "Segatella" results/biom_input/genus/*.kreport
grep -H "Segatella" results/bracken/genus/*.genus.bracken
grep -H "Segatella" results/bracken/species/*.species.bracken

# ---end--