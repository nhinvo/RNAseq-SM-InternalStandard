#!/usr/bin/env bash
#SBATCH --job-name=dge-smk
#SBATCH --time 10-0                         # wall clock limit
#SBATCH -p sched_mit_chisholm               # partition selection
#SBATCH -n 1                                # tasks (essentially threads)
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

source activate snakemake

snakemake -s workflow/Snakefile --unlock

# run job:
snakemake \
    --configfile "/nfs/chisholmlab001/kve/2022_chitin_bead_colonization_rnaseq/inputs/config.yaml" \
    -s workflow/Snakefile \
    --use-conda \
    --conda-frontend conda \
    --rerun-incomplete \
    --cluster "sbatch -p {resources.partition} -t {resources.time} -n {resources.ntasks} --mem={resources.mem} -o {resources.output} -e {resources.error}" \
    --jobs 10
    # --keep-going