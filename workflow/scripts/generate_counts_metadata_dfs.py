from pathlib import Path
import shutil
import pandas as pd
import json
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

# Snakemake
raw_gff_dir = Path(snakemake.input['raw_gff_dir'])
feature_types_to_keep = snakemake.config["feature_types_to_keep"]

log.info(f"snakemake.input['sample_counts']:\ntype:{type(snakemake.input['sample_counts'])}\ndata:{snakemake.input['sample_counts']}\n\n")

shutil.copy(snakemake.input['condition_table_path'], snakemake.output['samples'])
shutil.copy(snakemake.input["config_yaml_path"], snakemake.output["config"])

with open(snakemake.input['ref_json']) as file:
    ref_dict = json.load(file) 

# attributes and annotations

def generate_gff_df(gff_file):
    ## This method is heavily influenced from GFFpandas
    df = pd.read_csv(gff_file, sep="\t", comment="#",
            names=[
                "seq_id",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ],
    )
    attr_dict_series = df.attributes.apply(
        lambda attributes: dict(
            [
                key_value_pair.split(sep="=", maxsplit=1)
                for key_value_pair in attributes.split(";")
            ]
        )
    )
    key_set = set()
    attr_dict_series.apply(
        lambda at_dic: key_set.update(list(at_dic.keys()))
    )

    for attr in sorted(list(key_set)):
        df[attr] = attr_dict_series.apply(
            lambda attr_dict: attr_dict.get(attr)
        )

    return df
    

gffs = []
for gff_path in raw_gff_dir.iterdir():
    organism_name = str(gff_path.stem)
    attributes_df = generate_gff_df(gff_path)
    attributes_df['organism'] = organism_name
    gffs.append(attributes_df)

attributes_df = pd.concat(gffs)

log.info(f"attributes_df.columns:\n{attributes_df.columns}\n\n")

attributes_df = attributes_df.set_index('ID')

organism_type_ID_df = attributes_df[['organism', 'type']].copy()

attributes_df = attributes_df.set_index(['organism', 'type'], append=True)
attributes_df = attributes_df.reorder_levels(['organism', 'type', 'ID'])

raw_counts_df_list = []
raw_metadata_df_list = []
first_unique_ids = []

for i, path in enumerate(sorted(snakemake.input['sample_counts'])):
    path = Path(path)
    if path.suffix == '.tsv':
        # get sample ID from path
        sample_name = path.stem

        # read in HTseq TSV
        temp_df = pd.read_csv(path, sep='\t', names=['ID', sample_name])

        # check that IDs match
        if len(first_unique_ids) == 0:
            first_unique_ids = temp_df['ID'].unique()
        else:
            temp_unique_ids = temp_df['ID'].unique()
            assert first_unique_ids.all() == temp_unique_ids.all()

        temp_df = temp_df.set_index('ID')

        temp_metadata_df = temp_df[temp_df.index.str.contains('__')]

        temp_counts_df = temp_df[~temp_df.index.str.contains('__')]

        # append df to raw_counts_df_list
        raw_counts_df_list.append(temp_counts_df)
        raw_metadata_df_list.append(temp_metadata_df)

counts_df = pd.concat(raw_counts_df_list, axis=1)
metadata_df = pd.concat(raw_metadata_df_list, axis=1)

metadata_df.index = metadata_df.index.str.replace('__', '')

counts_df = counts_df.join(organism_type_ID_df)
counts_df = counts_df.set_index(['organism', 'type'], append=True)
counts_df = counts_df.reorder_levels(['organism', 'type', 'ID'])

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

counts_df['annotation_data', 'type'] = counts_df.index.get_level_values('type')
counts_df = counts_df.reset_index(level='type', drop=True)

counts_df.to_csv(Path(snakemake.output["counts"]), sep="\t")
metadata_df.to_csv(Path(snakemake.output["metadata"]), sep='\t')

annotation_data = counts_df['annotation_data']

annotation_dict = {}

for col, col_dict in ref_dict.items():
    annotation_dict[str(col)] = {}
    for i, term in col_dict.items():
        try:
            if i not in ['nan', "", 'None']:
                mi = annotation_data[annotation_data[col].str.contains(i, na=False)].index
                if term in ["",'None','nan']:
                    annotation_dict[str(col)][str(i)] = {org : [] for org in mi.get_level_values('organism').unique()}
                    for org, gene in mi.to_list():
                        annotation_dict[str(col)][str(i)][str(org)].append(gene)
                else:
                    annotation_dict[str(col)][str(term)] = {org : [] for org in mi.get_level_values('organism').unique()}
                    for org, gene in mi.to_list():
                        annotation_dict[str(col)][str(term)][str(org)].append(gene)
        except ValueError:
            log.debug(f"{col}\t{i}\t{term}")
            raise

cols2exclude = list(ref_dict.keys())
cols2exclude.extend(['start', 'end', 'attributes'])

for col in annotation_data.columns:
    if col not in cols2exclude:
        log.info(col)
        unique_values = annotation_data[col].unique()
        annotation_dict[str(col)] = {}
        for val in unique_values:
            mi = annotation_data[annotation_data[col]==val].index
            annotation_dict[str(col)][str(val)] = {org : [] for org in mi.get_level_values('organism').unique()}
            for org, gene in mi.to_list():
                annotation_dict[str(col)][str(val)][str(org)].append(gene)


with open(snakemake.output['annotation_json'], 'w') as outfile:
    json.dump(annotation_dict, outfile)

counts_data = counts_df['sample_data']
counts_dict = {}
for col in counts_data.columns:
    counts_dict[col] = {}
    for org in counts_data[col].index.get_level_values('organism').unique():
        counts_dict[col][org] = counts_data.loc[org][col].to_dict()

with open(snakemake.output['counts_json'], 'w') as outfile:
    json.dump(counts_dict, outfile)