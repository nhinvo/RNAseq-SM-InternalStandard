from pathlib import Path
import pandas as pd
import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

def create_attr_dict(attr_row):
    log.debug(attr_row)
    d = {}
    for key_value_pair in attr_row.split(";"):
        k_v_list = key_value_pair.split(sep="=", maxsplit=1)
        if len(k_v_list) == 2:
            k, v = k_v_list
            d[k] = v
    return d

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
        lambda attributes: create_attr_dict(attributes)
    )
    key_set = set()
    attr_dict_series.apply(lambda at_dic: key_set.update(list(at_dic.keys())))
    for attr in sorted(list(key_set)):
        df[attr] = attr_dict_series.apply(
            lambda attr_dict: attr_dict.get(attr)
        )
    return df

gffs = []
for genome, info in snakemake.config['input']['reference genomes'].items():
    attributes_df = generate_gff_df(info['annotation'])
    attributes_df['organism'] = genome
    gffs.append(attributes_df)
attributes_df = pd.concat(gffs)

# filtering by 'type' column 
feature_types_to_keep = snakemake.config["feature_types_to_count"]
attributes_df = attributes_df[attributes_df['type'].isin(feature_types_to_keep)]

attributes_df = attributes_df.set_index(['organism','ID'])
attributes_df.to_csv(snakemake.output[0], sep='\t')
