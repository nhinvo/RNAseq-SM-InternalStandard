from pathlib import Path
import pandas as pd
import logging as log
import urllib.error
import Bio.KEGG.REST as kegg_rest
import json

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

raw_gff_dir = Path(snakemake.input['raw_gff_dir']) # Path('/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data/culture_genome_annotations')

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

attributes_df = attributes_df.rename(columns={'ID':'long_ID'})

log.info(f"attributes_df.columns:\n{attributes_df.columns}\n\n")

attributes_df = attributes_df.set_index(['organism', 'type', 'long_ID'])

kegg_pathway_columns = {'em_BRITE':'brite', 'em_KEGG_Module':'module', 'em_KEGG_Pathway':'pathway', 'em_KEGG_Reaction':'reaction', 'em_KEGG_rclass':'rclass'}

def ident_kegg_term(db, id):
    try: 
        line = kegg_rest.kegg_find(db, id).readline()
        return line.split('\n')[0].split('\t')[1]
    except IndexError as e:
        log.error(f'problematic term:\n{kegg_rest.kegg_find(db, id).readline()}')
        return ""
    except urllib.error.URLError as e:
        return e.read().decode("utf8", 'ignore')

ref_dict = {}
for col, db in kegg_pathway_columns.items():
    log.info(f"{col}\t{db}")
    unique_terms = set()
    ref_dict[col] = {}
    if col in attributes_df.columns:
        for row in attributes_df[col].to_list():
            unique_terms |= set(str(row).split(","))
        unique_terms -= {'-', 'nan'}
        term_table = []
        for i in unique_terms:
            term = ident_kegg_term(db, i)
            log.info(f"{i}\t{term}")
            ref_dict[col][i] = term

with open(snakemake.output['ref_json'], "w") as outfile:
# with open('/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/results_bwa/ref.json', "w") as outfile:
    json.dump(ref_dict, outfile)