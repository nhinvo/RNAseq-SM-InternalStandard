from pathlib import Path
import pandas as pd
import numpy as np
import json
import yaml

# rpy2 imports
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter

import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

base = rpackages.importr('base')
utils = rpackages.importr('utils')
deseq2 = rpackages.importr('DESeq2')

def get_dge_table(results_table):

    def symlog10(x):
        if x > 0:
            return np.log10(x+1)
        elif x < 0:
            return -np.log10(-x+1)
        else:
            return 0
        
    results_table['symlog10baseMean'] = results_table['baseMean'].apply(lambda x: symlog10(x))

    cols = list(results_table.columns.values)
    
    index_of_baseMean = cols.index('baseMean')
    cols.remove('symlog10baseMean')
    cols.insert(index_of_baseMean+1, 'symlog10baseMean')

    results_table = results_table[cols]

    return results_table

def main():

    counts_df = pd.read_csv(Path(snakemake.input['counts']), sep='\t', index_col=[0,1], header=[0,1])

    with open(Path(snakemake.input['comparison_yaml'])) as infile:
        comparisons = yaml.load(infile, Loader=yaml.FullLoader)

    log.debug(f"counts_df:\n{counts_df}")

    counts_df = counts_df.reset_index(level=1, col_fill='annotation_data')

    log.debug(f"counts_df after reindexing:\n{counts_df}\n\n")

    counts_df = counts_df['sample_data']

    for organism in comparisons.keys():
        for comparison in comparisons[organism].keys():

            log.info(f"organism: {organism}\t comparison: {comparison}")

            sample_use = pd.Series(comparisons[organism][comparison]['sample_use'], name='sample_use')
            sample_use = sample_use[sample_use != 'none']

            conditions = list(sample_use.unique())

            return_results_df = len(conditions) == 2 and 'control' in conditions and 'treatment' in conditions

            log.info(f'{organism}\n')

            included_samples = list(sample_use.index.values)

            comp_count_df = counts_df.loc[organism][included_samples]

            sample_use_df = pd.DataFrame(sample_use)

            log.debug(f"comp_count_df:\n{comp_count_df}")
            log.debug(f"sample_use_df:\n{sample_use_df}")

            # organism must be present in all samples. Good test to tell if it is from experience...
            if min(comp_count_df.mean()) > 1:

                comparisons[organism][comparison]['status'] = 'running'

                # DEseq process
                # 1a. transfer count df into r df
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    r_comp_count_df = robjects.conversion.py2rpy(comp_count_df)

                robjects.globalenv['r_comp_count_df'] = r_comp_count_df
                
                # 1b. transfer condition df into r df
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    r_comparison_condition_df = robjects.conversion.py2rpy(sample_use_df)
                robjects.globalenv['r_comparison_condition_df'] = r_comparison_condition_df
                
                # 2. create DESeqDataSet object
                # exp design
                robjects.r(f"""dds <- DESeqDataSetFromMatrix(countData = r_comp_count_df, 
                                                            colData = r_comparison_condition_df, 
                                                            design = ~ sample_use)""")
                
                # 3. run DEseq command
                dds_processed = robjects.r("DESeq(dds)")
                robjects.globalenv['dds_processed'] = dds_processed

                # get normalized counts df
                r_normalized_counts_df = robjects.r("counts(dds_processed, normalized=TRUE)")
                
                if return_results_df:
                    # 4a. set up comparison controls and treatments
                    contrast_string_list = robjects.StrVector(['sample_use', 'control', 'treatment'])

                    # 4b. get results df
                    r_results = deseq2.results(dds_processed, contrast = contrast_string_list)
                    robjects.globalenv['r_results'] = r_results
                    r_results_df = robjects.r("as.data.frame(r_results)")

                # 5. get rlog and vsd dfs
                rlog_output = deseq2.rlog(dds_processed, blind=False)
                robjects.globalenv['rlog_output'] = rlog_output
                r_rlog_df = robjects.r("assay(rlog_output)")

                robjects.r("vst_output <- varianceStabilizingTransformation(dds, blind=FALSE)")
                r_vst_df = robjects.r("assay(vst_output)")
                
                # 6. transfer normalized counts, rlog, and vst df to pandas
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    rlog_array = robjects.conversion.rpy2py(r_rlog_df)
                
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    vst_array = robjects.conversion.rpy2py(r_vst_df)
                
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    normalized_counts_array = robjects.conversion.rpy2py(r_normalized_counts_df)

                if return_results_df:
                    # 7. transfer results df to pandas
                    with localconverter(robjects.default_converter + pandas2ri.converter):
                        results_table = robjects.conversion.rpy2py(r_results_df)
                        
                    results_table = results_table.rename({'seq_id':'ID'})
                    results_table['ID'] = comp_count_df.index.values
                    results_table = results_table.set_index('ID')

                normalized_counts_df = pd.DataFrame(normalized_counts_array, index=comp_count_df.index, columns=included_samples)
                rlog_df = pd.DataFrame(rlog_array, index=comp_count_df.index, columns=included_samples)
                vst_df = pd.DataFrame(vst_array, index=comp_count_df.index, columns=included_samples)

                comparisons[organism][comparison]['rlog'] = rlog_df.to_dict(orient='tight')
                
                # post-DEseq2 analysis
                if return_results_df:
                    all_dge_table = get_dge_table(results_table)
                    log.debug(f'{all_dge_table}\n')
                    comparisons[organism][comparison]['results'] = all_dge_table.to_dict()
                
                comparisons[organism][comparison]['status'] = 'finished'

            else:
                comparisons[organism][comparison]['status'] = 'failed QC'

    with open(Path(snakemake.output['comparison_data']), 'w') as outfile:
        json.dump(comparisons, outfile)
                
if __name__ == '__main__':
    main()