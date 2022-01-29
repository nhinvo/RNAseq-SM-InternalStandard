from pathlib import Path
import shutil
from matplotlib.colors import same_color
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd

import logging as log

# Snakemake
feature_types_to_keep = snakemake.config["results"]["feature_types_to_keep"]
results_dir = Path(snakemake.config["results"])
htseq2dir = Path(snakemake.config["output"]["feature_count"])
raw_gff_dir = Path(snakemake.config["input"]["gff_refs"])
condition_table_path = Path(snakemake.config["samples"])

# Non-Snakemake
# with open('../../config/config.yaml', 'r') as file:
#     config = yaml.safe_load(file)

# results_dir = Path(config["results"])
# htseq2dir = Path(config["output"]["feature_count"])
# raw_gff_dir = Path(config["input"]["gff_refs"])
# condition_table_path = Path('../../config/samples.tsv')

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

def makeOutDir(outputDir, folderName):
    outdir = outputDir / folderName
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir

comparisons_dir = makeOutDir(results_dir, 'comparisons')
metadata_figures_dir = makeOutDir(results_dir, 'metadata_figures')
metadata_tables_dir = makeOutDir(results_dir, 'metadata_tables')

app_dfs_dir = makeOutDir(results_dir, 'app_dfs')
shutil.copy(condition_table_path, app_dfs_dir / condition_table_path.name)

# attributes and annotations

gffs = []
for gff_path in gff_dir.iterdir():
    organism_name = str(gff_path.stem)
    annotation = gffpd.read_gff3(gff_path)
    attributes_df = annotation.attributes_to_columns()
    attributes_df['organism'] = organism_name
    gffs.append(attributes_df)

attributes_df = pd.concat(gffs)

attributes_df = attributes_df.rename(columns={'ID':'long_ID'})

attributes_df = attributes_df.set_index('long_ID')

organism_type_ID_df = attributes_df[['organism', 'type']]

organism_type_ID_df.to_csv(results_dir / "organism_type_ID_df.tsv", sep="\t")

raw_counts_df_list = []
raw_metadata_df_list = []
first_unique_ids = []

for i, path in enumerate(sorted(list(htseq2dir.iterdir()))):
    if path.suffix == '.tsv':
        # get sample ID from path
        sample_name = path.stem

        # read in HTseq TSV
        temp_df = pd.read_csv(path, sep='\t', names=['long_ID', sample_name])

        # check that long_IDs match
        if len(first_unique_ids) == 0:
            first_unique_ids = temp_df['long_ID'].unique()
        else:
            temp_unique_ids = temp_df['long_ID'].unique()
            assert first_unique_ids.all() == temp_unique_ids.all()

        temp_df = temp_df.set_index('long_ID')

        temp_metadata_df = temp_df[temp_df.index.str.contains('__')]

        temp_counts_df = temp_df[~temp_df.index.str.contains('__')]

        # append df to raw_counts_df_list
        raw_counts_df_list.append(temp_counts_df)
        raw_metadata_df_list.append(temp_metadata_df)

counts_df = pd.concat(raw_counts_df_list, axis=1)
counts_df = counts_df.add_prefix('sample_')

metadata_df = pd.concat(raw_metadata_df_list, axis=1)
metadata_df = metadata_df.add_prefix('sample_')
metadata_df.index = metadata_df.index.str.replace('__', '')

counts_df = counts_df.join(organism_type_ID_df)
counts_df = counts_df.set_index(['organism', 'type'], append=True)
counts_df = counts_df.reorder_levels(['organism', 'type', 'long_ID'])

if feature_types_to_keep:
    counts_df = counts_df[counts_df.index.get_level_values('type').isin(feature_types_to_keep)]

feature_df = counts_df.groupby(['type']).sum() 
metadata_df = pd.concat([feature_df, metadata_df])


# todo: combine counts and attributes dfs 
counts_df.to_csv(results_dir / "counts.tsv", sep="\t")
attributes_df.to_csv(results_dir / "attributes.tsv", sep="\t")
metadata_df.to_csv(results_dir / 'metadata.tsv', sep='\t')


