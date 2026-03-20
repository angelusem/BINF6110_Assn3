#!/bin/bash
#SBATCH --job-name=a3_fasterq
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/%u/assn3_binf6110/logs/fasterq_%j.out
#SBATCH --error=/scratch/%u/assn3_binf6110/logs/fasterq_%j.err

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
