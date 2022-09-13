from pathlib import Path
import pandas as pd
import logging as log
import json
import yaml

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

metadata_df = pd.read_csv(snakemake.input['metadata'], sep='\t', header=0, index_col=0)
samples_df = pd.read_csv(snakemake.input['samples'], sep='\t', header=0)

log.debug(f"metadata_df:\n{metadata_df}\n\n")
log.debug(f"samples_df:\n{samples_df}\n\n")

all_data = {}

with open(snakemake.input['config']) as file:
    all_data['config'] = yaml.load(file, Loader=yaml.FullLoader)

with open(snakemake.input['comparison_data']) as file:
    all_data['comparisons'] = json.load(file) 

with open(snakemake.input['counts_json']) as file:
    all_data['count_data'] = json.load(file) 

with open(snakemake.input['annotation_json']) as file:
    all_data['annotation_data'] = json.load(file) 

all_data['samples'] = samples_df.to_dict(orient='records')

all_data['metadata'] = metadata_df.to_dict(orient='dict')

with open(snakemake.output['samples'], "w") as outfile:
    json.dump(all_data['samples'], outfile)

with open(snakemake.output['metadata'], "w") as outfile:
    json.dump(all_data['metadata'], outfile)

with open(snakemake.output['config'], "w") as outfile:
    json.dump(all_data['config'], outfile)

with open(snakemake.output['data_json'], "w") as outfile:
    json.dump(all_data, outfile)

