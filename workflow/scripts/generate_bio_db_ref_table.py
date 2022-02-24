from pathlib import Path
import pandas as pd
import logging as log
import urllib.error
import gffpandas.gffpandas as gffpd
import Bio.KEGG.REST as kegg_rest

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

raw_gff_dir = Path("/nfs/chisholmlab001/kve/2021_Sar11Pro_RNAseq_Project/data/input_data/culture_genome_annotations") #snakemake.input['raw_gff_dir'])

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

attributes_df = attributes_df.set_index(['organism', 'type', 'long_ID'])

kegg_pathway_columns = {'em_BRITE':'brite', 'em_KEGG_Module':'module', 'em_KEGG_Pathway':'pathway', 'em_KEGG_Reaction':'reaction', 'em_KEGG_rclass':'rclass'}

def ident_kegg_term(db, id):
    try: 
        return kegg_rest.kegg_find(db, id).readline().split('\n')[0].split('\t')[1]
    except IndexError as e:
        log.error(f'problematic term:\n{kegg_rest.kegg_find(db, id).readline()}')
        return ""
    except urllib.error.URLError as e:
        return e.read().decode("utf8", 'ignore')

dfs = []
for col, db in kegg_pathway_columns.items():
    log.info(f"{col}\t{db}")
    unique_terms = set()
    if col in attributes_df.columns:
        for row in attributes_df[col].to_list():
            unique_terms |= set(str(row).split(","))
        unique_terms -= {'-', 'nan'}
        term_table = []
        for i in unique_terms:
            term = ident_kegg_term(db, i)
            log.info(f"{col}\t{i}\t{term}")
            term_table.append([col, i, term])
        dfs.append(pd.DataFrame(term_table, columns=['ref_col','id','term']))

ref_table = pd.concat(dfs, axis=0)
ref_table.to_csv(snakemake.output['kegg_refs'], sep='\t', index=False)
