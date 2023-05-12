from pathlib import Path
import pandas as pd
import logging as log
import json
import yaml
import gzip

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

open(snakemake.output[0],'w').close()

mapping_metadata = pd.read_table(snakemake.input['mapping_metadata'], index_col=0)
organism_occurance = pd.read_table(snakemake.input['organism_occurance'], index_col='organism')
gene_sparsity = pd.read_table(snakemake.input['gene_sparsity'], index_col='organism')
samples = pd.read_table(snakemake.config["input"]["sample table"], index_col="sample")
annotations = pd.read_table(snakemake.input['annotations'], index_col=['organism', 'ID'])
raw_counts = pd.read_table(snakemake.input['counts'], index_col=['organism', 'ID'])
bio_ref_df = pd.read_table(snakemake.input['bio_df_ref'])

data = {
    'metadata': {
        'mapping quality': mapping_metadata.to_dict(),
        'organism occurance': organism_occurance.to_dict(),
        'gene sparsity': gene_sparsity.to_dict(),
    },
    'samples' : samples.to_dict(),
    'config' : snakemake.config,
    'annotations': {
        org: annotations.loc[org].to_dict()
        for org in annotations.index.get_level_values('organism').unique()
    },
    'raw_counts': {
        org: raw_counts.loc[org].to_dict()
        for org in raw_counts.index.get_level_values('organism').unique()
    },
    'id ref table': bio_ref_df.to_dict(),
    'deseq2' : {}
}

for rlog, vst, result in zip(snakemake.input['rlogs'], snakemake.input['vsts'], snakemake.input['deseq2_results']):
    assert Path(rlog).parent == Path(vst).parent == Path(result).parent
    org = Path(rlog).parent.parent.parent.name
    exp = Path(rlog).parent.parent.name
    comp = Path(rlog).parent.name
    if org not in data['deseq2'].keys():
        data['deseq2'][org] = {}
    if exp not in data['deseq2'][org].keys():
        data['deseq2'][org][exp] = {}
    data['deseq2'][org][exp][comp] = {
        'rlog' : pd.read_table(rlog, index_col='ID').to_dict(),
        'vst' : pd.read_table(vst, index_col='ID').to_dict(),
        'results' : pd.read_table(result, index_col='ID').to_dict()
    }

with gzip.open(snakemake.output[0], 'w') as fout:
    fout.write(json.dumps(data).encode('utf-8'))  