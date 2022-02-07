from pathlib import Path
import shutil
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
import math

# rpy2 imports
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter

# R package names
packnames = ('base', 'utils', 'DESeq2')

# R vector of strings
from rpy2.robjects.vectors import StrVector

# Selectively install what needs to be install.
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

base = rpackages.importr('base')
utils = rpackages.importr('utils')
deseq2 = rpackages.importr('DESeq2')

import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.DEBUG)

def get_dge_table(results_table):

    def symlog10(x):
        if x > 0:
            return math.log10(x+1)
        elif x < 0:
            return -math.log10(-x+1)
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

    counts_df = pd.read_csv(Path(snakemake.input['counts']), sep='\t', index_col=[0,1,2], header=[0,1])
    conditions_df =  pd.read_csv(Path(snakemake.input['samples']), sep='\t', index_col='sample_name')

    log.debug(f"counts_df:\n{counts_df}")

    # extracts columns that include the word "comparison"
    comparisons_include_df = conditions_df[conditions_df.columns[np.array(conditions_df.columns.str.contains('comparison'))]]

    # creates new dataframe from comparisons_include_df that replaces control/treatment with True in order to subset future datasets
    conditions_df_sans_comparisons = conditions_df[conditions_df.columns[~np.array(conditions_df.columns.str.contains('comparison'))]]

    counts_df = counts_df.reset_index(level=1, col_fill='annotation_data')

    log.debug(f"counts_df after reindexing:\n{counts_df}\n\n")

    annotation_data = counts_df[['annotation_data']]

    log.debug(f"annotation_data:\n{annotation_data}")

    counts_df = counts_df['sample_data']

    rlog_dfs = []
    results_dfs = [annotation_data]

    for comparison in comparisons_include_df.columns:
        
        log.info(comparison)

        included_samples_series = comparisons_include_df[comparisons_include_df[comparison].notna()][comparison]

        conditions = list(included_samples_series.unique())
        
        return_results_df = len(conditions) == 2 and 'control' in conditions and 'treatment' in conditions

        log.info(f"conditions {return_results_df} for {conditions}")

        comparison_condition_df = conditions_df_sans_comparisons.loc[included_samples_series.index]
        comparison_condition_df[comparison] = included_samples_series

        comparison_samples = included_samples_series.index.values 

        comparison_counts_df = counts_df[[f"sample_{s:02d}" for s in comparison_samples]]

        organisms_result_dfs = []
        for organism in comparison_counts_df.index.unique(0):
            log.info(f'{organism}\n')

            comparison_organism_count_df = comparison_counts_df.loc[organism, :]

            # organism must be present in all samples. Good test to tell if it is from experience...
            # log.info('{} - organism lowest sample read mapping in comparison: {}'.format(organism, min(comparison_organism_count_df.mean())))
            if min(comparison_organism_count_df.mean()) > 1:

                # DEseq process
                # 1a. transfer count df into r df
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    r_comparison_organism_count_df = robjects.conversion.py2rpy(comparison_organism_count_df)

                robjects.globalenv['r_comparison_organism_count_df'] = r_comparison_organism_count_df
                
                # 1b. transfer condition df into r df
                with localconverter(robjects.default_converter + pandas2ri.converter):
                    r_comparison_condition_df = robjects.conversion.py2rpy(comparison_condition_df)
                robjects.globalenv['r_comparison_condition_df'] = r_comparison_condition_df
                
                # 2. create DESeqDataSet object
                # exp design
                robjects.r(f"""dds <- DESeqDataSetFromMatrix(countData = r_comparison_organism_count_df, 
                                                            colData = r_comparison_condition_df, 
                                                            design = ~ {comparison})""")
                
                # 3. run DEseq command
                dds_processed = robjects.r("DESeq(dds)")
                robjects.globalenv['dds_processed'] = dds_processed

                # get normalized counts df
                r_normalized_counts_df = robjects.r("counts(dds_processed, normalized=TRUE)")
                
                if return_results_df:
                    # 4a. set up comparison controls and treatments
                    contrast_string_list = robjects.StrVector([comparison, 'control', 'treatment'])

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
                        
                    results_table = results_table.rename({'seq_id':'long_ID'})
                    results_table['long_ID'] = comparison_organism_count_df.index.values
                    results_table = results_table.set_index('long_ID')

                normalized_counts_df = pd.DataFrame(normalized_counts_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                rlog_df = pd.DataFrame(rlog_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                vst_df = pd.DataFrame(vst_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                

                # post-DEseq2 analysis

                if return_results_df:

                    all_dge_table = get_dge_table(results_table)
                    
                    #results df
                    multiindex_results_df = results_table
                    multiindex_results_df.columns = pd.MultiIndex.from_product([[comparison], results_table.columns])
                    multiindex_results_df.index = pd.MultiIndex.from_product([[organism], results_table.index])
                    # log.debug(f"multiindex_results_df.shape before filter:\n{multiindex_results_df.shape}")
                    # multiindex_results_df = multiindex_results_df.loc[~multiindex_results_df.index.duplicated(keep='first')]
                    # log.debug(f"multiindex_results_df.shape after filter:\n{multiindex_results_df.shape}")
                    log.debug(f"multiindex_results_df:\n{multiindex_results_df}")
                    organisms_result_dfs.append(multiindex_results_df)

                    #rlog df
                    multiindex_rlog_df = rlog_df
                    multiindex_rlog_df.columns = pd.MultiIndex.from_product([[organism], [comparison], rlog_df.columns])
                    multiindex_rlog_df.index = pd.MultiIndex.from_product([[organism], rlog_df.index])
                    log.debug(f"multiindex_rlog_df:\n{multiindex_rlog_df}")
                    rlog_dfs.append(multiindex_rlog_df)
        
        if len(organisms_result_dfs) == 1:
            results_dfs.append(organisms_result_dfs[0])
        if len(organisms_result_dfs) > 1:
            results_df = pd.concat(organisms_result_dfs, axis=0)
            results_dfs.append(results_df)
        

    app_results_df = pd.concat(results_dfs, axis=1)
    app_rlog_df = pd.concat(rlog_dfs, axis=1)

    log.debug(f"app_results_df:\n{app_results_df}")
    log.debug(f"app_rlog_df:\n{app_rlog_df}")

    app_results_df.to_csv(Path(snakemake.output['results']), sep='\t')
    app_rlog_df.to_csv(Path(snakemake.output['rlog']), sep='\t')


if __name__ == '__main__':
    main()