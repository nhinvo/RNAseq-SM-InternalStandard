from pathlib import Path
import pandas as pd
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

attributes_df = pd.read_table(snakemake.input['annotations'], index_col=[0,1])

raw_counts_df_list = []
raw_metadata_df_list = []
reference_index = None

for path in [Path(p) for p in snakemake.input['counts']]:
    sample_name = path.stem
    temp_df = pd.read_table(path, names=['ID', sample_name])
    temp_df = temp_df.set_index('ID')

    if isinstance(reference_index, pd.Index):
        assert reference_index.all() == temp_df.index.all()
    else:
        reference_index = temp_df.index

    temp_metadata_df = temp_df[temp_df.index.str.contains('__')]
    temp_counts_df = temp_df[~temp_df.index.str.contains('__')]

    # append df to raw_counts_df_list
    raw_counts_df_list.append(temp_counts_df)
    raw_metadata_df_list.append(temp_metadata_df)

counts_df = pd.concat(raw_counts_df_list, axis=1)
metadata_df = pd.concat(raw_metadata_df_list, axis=1)

metadata_df.index = metadata_df.index.str.replace('__', '')

print(counts_df)
print(attributes_df)
counts_df = counts_df.loc[attributes_df.index.get_level_values('ID')]
counts_df.index = attributes_df.index

feature_df = counts_df.join(attributes_df['type']).groupby(['type']).sum() 
metadata_df = pd.concat([feature_df, metadata_df])

counts_df.to_csv(Path(snakemake.output["raw_counts"]), sep="\t")
metadata_df.to_csv(Path(snakemake.output["mapping_metadata"]), sep='\t')
counts_df.join(attributes_df).to_csv(Path(snakemake.output["raw_counts_w_annotations"]), sep="\t")

organism_occurance = counts_df.groupby(level='organism').sum().div(counts_df.sum())
organism_occurance.to_csv(Path(snakemake.output["organism_occurance"]), sep='\t')

sparsity_df = (counts_df >= 1).groupby(level='organism').sum().div(counts_df.index.get_level_values('organism').value_counts(), axis='index')
sparsity_df.index.name = 'organism'
sparsity_df.to_csv(Path(snakemake.output["gene_sparsity"]), sep="\t")

