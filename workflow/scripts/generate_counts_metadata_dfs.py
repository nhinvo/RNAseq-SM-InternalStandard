from pathlib import Path
import shutil
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd

import logging as log

# Snakemake

htseq2dir = Path(snakemake.input['sample_counts'][0]).parent
raw_gff_dir = Path(snakemake.input['raw_gff_dir'])
feature_types_to_keep = snakemake.config["feature_types_to_keep"]

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

shutil.copy(snakemake.input['condition_table_path'], snakemake.output['samples'])
shutil.copy(snakemake.input["config_yaml_path"], snakemake.output["config"])

# attributes and annotations

gffs = []
for gff_path in raw_gff_dir.iterdir():
    organism_name = str(gff_path.stem)
    annotation = gffpd.read_gff3(gff_path)
    attributes_df = annotation.attributes_to_columns()
    attributes_df['organism'] = organism_name
    gffs.append(attributes_df)

attributes_df = pd.concat(gffs)

attributes_df = attributes_df.rename(columns={'ID':'long_ID'})

log.info(f"attributes_df.columns:\n{attributes_df.columns}\n\n")

attributes_df = attributes_df.set_index('long_ID')

organism_type_ID_df = attributes_df[['organism', 'type']].copy()

attributes_df = attributes_df.set_index(['organism', 'type'], append=True)
attributes_df = attributes_df.reorder_levels(['organism', 'type', 'long_ID'])

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

counts_df.columns = pd.MultiIndex.from_product([['sample_data'], counts_df.columns])
attributes_df.columns = pd.MultiIndex.from_product([['annotation_data'], attributes_df.columns])
counts_df = counts_df.join(attributes_df) 

log.info(f"counts_df:\n{counts_df}\n\n")
log.info(f"metadata_df:\n{metadata_df}\n\n")

counts_df = counts_df.loc[~counts_df.index.duplicated(keep='first')]

counts_df.to_csv(Path(snakemake.output["counts"]), sep="\t")
metadata_df.to_csv(Path(snakemake.output["metadata"]), sep='\t')
