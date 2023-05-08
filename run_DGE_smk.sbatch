#!/usr/bin/env bash
#SBATCH --job-name=dge-smk
#SBATCH --time 5-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

source activate snakemake

# run job:
snakemake --profile profile
