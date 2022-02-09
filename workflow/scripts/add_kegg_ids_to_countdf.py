from pathlib import Path
import pandas as pd
import logging as log
import urllib.error

import Bio.KEGG.REST as kegg_rest

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

counts_df = pd.read_csv(Path(snakemake.input['counts']), sep='\t', index_col=[0,1,2], header=[0,1])

log.debug(f"counts_df:\n{counts_df}\n\n")

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
    if col in counts_df['annotation_data'].columns:
        for row in counts_df['annotation_data'][col].to_list():
            unique_terms |= set(str(row).split(","))
        unique_terms -= {'-', 'nan'}
        term_table = []
        for i in unique_terms:
            term = ident_kegg_term(db, i)
            log.info(f"{col}\t{i}\t{term}")
            term_table.append([col, db, i, term])
        dfs.append(pd.DataFrame(term_table, columns=['ref_col', 'kegg_df_ref', 'kegg_id', 'kegg_term']))

ref_table = pd.concat(dfs, axis=0)
ref_table.to_csv(snakemake.output['kegg_refs'], sep='\t', index=False)

counts_df_kegg_cols = ref_table.ref_col.unique()

counts_df['annotation_data'] = counts_df['annotation_data'].fillna(value={col:"" for col in counts_df_kegg_cols})
counts_df['annotation_data'] = counts_df['annotation_data'].replace(to_replace="-", value={col:"" for col in counts_df_kegg_cols})

for col in ref_table.ref_col.unique():
    ref_table_subset = ref_table[ref_table['ref_col']==col]
    ref_dict = {row.loc['kegg_id']:row.loc['kegg_term'] for i, row in ref_table_subset.iterrows()}
    # log.debug(f"ref_dict:\n{ref_dict}")
    col_ids = counts_df['annotation_data'][col].to_list()
    col_terms = []
    for row in col_ids:
        if row == "":
            row_label = ""
        else:
            row_label = ", ".join([str(ref_dict[i]) for i in str(row).split(',')])
        col_terms.append(row_label)

    counts_df.loc[:,('annotation_data',f"{col}_terms")] = col_terms

counts_df.to_csv(Path(snakemake.input['counts']), sep='\t')
