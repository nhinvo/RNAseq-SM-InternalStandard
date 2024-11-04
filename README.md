# RNAseq Pipeline with Internal Standard Correction 
## Installations

## Snakemake Setup 
### Create samples.tsv file: 
Required columns:
1. "sample" - Unique sample name 
2. "forward read" - Absolute path to forward reads .fastq file 
3. "reverse read" - Absolute path to reverse reads .fastq file

### Create internal_standard_concentration.tsv file:
Required columns:
1. "standard_name" - Unique name of internal standard added 
2. "standard_group" - Group that standard is in (note: future edits to remove this)
3. "concentration (ng/ul)" - Concentration of standard added 
4. "volume_added (ul)" - Volume of standard added 

### Create cell_count.tsv file:
Required columns:
1. "sample" - Unique sample name. Should match with "sample" column from samples.tsv	
2. "total_cell_count" - Count of cells in sample 

### Edit config.yaml file:
1. Edit names (yaml keys) of reference genomes and their paths
2. Edit path to folder to store intermediate files: "scratch directory"
3. Edit length of sequenced read (from fastq file): "read length"
4. Edit minimum number of samples a standard has to be in: "minimum standard sample"

### Edit profile/config.yaml file:
1. Edit partition name in "default-resources" - "partition"
2. Edit any other resources as needed 

### Edit run_RNAseq_SM.sbatch file: 
1. Edit slurm SBATCH resource specifications as needed 

## Running RNAseq Snamake Pipeline
1. Run pipeline by: `sbatch run_RNAseq_SM.sbatch`
    - Note: create logs/ folder before submitting job 
2. Check log files in logs/ folder 