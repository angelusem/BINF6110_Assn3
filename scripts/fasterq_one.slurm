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
