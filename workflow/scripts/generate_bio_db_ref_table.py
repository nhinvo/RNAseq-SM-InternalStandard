from pathlib import Path
import pandas as pd
import logging as log
import urllib.error
import Bio.KEGG.REST as kegg_rest

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

kegg_pathway_columns = {'em_BRITE':'brite', 'em_KEGG_Module':'module', 'em_KEGG_Pathway':'pathway', 'em_KEGG_Reaction':'reaction', 'em_KEGG_rclass':'rclass'}

kegg_ids = {k:set() for k in kegg_pathway_columns.values()}

for genome, info in snakemake.config['input']['reference genomes'].items():
    attributes = pd.read_table(info['annotation'], header=None, comment="#")[8].to_list()
    for gene in attributes:
        for annot in gene.split(";"):
            k, v = annot.split("=", maxsplit=1)
            if k in kegg_pathway_columns.keys():
                kegg_ids[kegg_pathway_columns[k]].update(v.split(","))

def ident_kegg_term(db, id):
    try: 
        line = kegg_rest.kegg_find(db, id).readline()
        return line.split('\n')[0].split('\t')[1]
    except IndexError as e:
        log.error(f'problematic term {id}:\n{kegg_rest.kegg_find(db, id).readline()}')
        return ""
    except urllib.error.URLError as e:
        return e.read().decode("utf8", 'ignore')

ref_table = []
for col, db in kegg_pathway_columns.items():
    s = kegg_ids[db]
    s -= {'-', 'nan'}
    log.info(f"{col} with {len(s)} entries")
    for i in s:
        log.info(i)
        ref_table.append([col, db, i, ident_kegg_term(db, i)])

pd.DataFrame(ref_table, columns=['eggnog column', 'kegg hierarchy id', 'kegg id', 'kegg name']).to_csv(snakemake.output[0], sep='\t', index=False)
