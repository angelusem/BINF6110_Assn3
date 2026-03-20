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
