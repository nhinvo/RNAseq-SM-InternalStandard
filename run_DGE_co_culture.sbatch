#!/usr/bin/env bash
# 
#SBATCH --job-name=DGE-co-cul
#
# Number of tasks/cores for job
#SBATCH -n 10
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 72:00:00
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
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
#SBATCH --mem=250G
#
#SBATCH -o logs/%j_slurm_output.txt
#SBATCH -e logs/%j_slurm_error.txt

source activate DGE_smk_minimal

SNAKEFILE=workflow/Snakefile

# lints and formats snakefile
snakemake --snakefile $SNAKEFILE --lint 
snakemake --snakefile $SNAKEFILE --unlock
snakefmt $SNAKEFILE

# examine full workflow using:
snakemake --forceall -n --snakefile $SNAKEFILE --dag > all_rules_dag.dot

# run job:
snakemake --cores all --resources mem_mb=250000 --snakefile $SNAKEFILE --use-conda --conda-frontend conda --rerun-incomplete
# get report once job is finished
snakemake --snakefile $SNAKEFILE --report report.html