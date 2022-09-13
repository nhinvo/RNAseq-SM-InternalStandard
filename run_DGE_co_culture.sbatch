#!/usr/bin/env bash
# 
#SBATCH --job-name=DGE-co-cul
#
# Number of tasks/cores for job
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --mem 250000
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 5-0
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
#SBATCH --comment="Diff expression analysis of co-culture"
#
# emails all notifications
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kve@mit.edu
#
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate snakemake

SNAKEFILE=workflow/Snakefile

# lints and formats snakefile
snakemake --snakefile $SNAKEFILE --lint 
snakemake --snakefile $SNAKEFILE --unlock
snakefmt $SNAKEFILE

# run job:
snakemake --cores all --resources mem_mb=500000 --snakefile $SNAKEFILE --use-conda --conda-frontend conda --rerun-incomplete
