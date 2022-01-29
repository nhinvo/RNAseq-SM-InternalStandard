from pathlib import Path
import shutil
from matplotlib.colors import same_color
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()
from sklearn.decomposition import PCA
import gffpandas.gffpandas as gffpd
import math

# rpy2 imports
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.conversion import localconverter
base = importr('base')
utils = importr('utils')
deseq2 = importr('DESeq2')

import logging as log

log.basicConfig(format='%(levelname)s:%(message)s', level=log.INFO)

def makeOutDir(outputDir, folderName):
    """
    Makes an out directory if there is not one. Returns file path as Path object
    Inputs:
    outputDir - Pathlib directory object that directory will be created within
    folderName - Name of directory to be made/returned
    """
    outdir = outputDir / folderName
    if not outdir.exists():
        outdir.mkdir()
    return outdir

def get_htseq2_and_metadata_df(htseq2dir, organism_type_ID_df, feature_types_to_keep=None):
    raw_counts_df_list = []
    raw_metadata_df_list = []
    first_unique_ids = []

    for i, path in enumerate(sorted(list(htseq2dir.iterdir()))):
        if path.suffix == '.tsv':
            # get sample ID from path
            sample_name = path.stem

            # read in HTseq TSV
            temp_df = pd.read_csv(path, sep='\t', names=['long_ID', sample_name])

            # check that long_IDs match
            if len(first_unique_ids) == 0:
                first_unique_ids = temp_df['long_ID'].unique()
            else:
                temp_unique_ids = temp_df['long_ID'].unique()
                assert first_unique_ids.all() == temp_unique_ids.all()

            temp_df = temp_df.set_index('long_ID')

            temp_metadata_df = temp_df[temp_df.index.str.contains('__')]

            temp_counts_df = temp_df[~temp_df.index.str.contains('__')]

            # append df to raw_counts_df_list
            raw_counts_df_list.append(temp_counts_df)
            raw_metadata_df_list.append(temp_metadata_df)

    counts_df = pd.concat(raw_counts_df_list, axis=1)
    counts_df = counts_df.add_prefix('sample_')

    metadata_df = pd.concat(raw_metadata_df_list, axis=1)
    metadata_df = metadata_df.add_prefix('sample_')
    metadata_df.index = metadata_df.index.str.replace('__', '')

    counts_df = counts_df.join(organism_type_ID_df)
    counts_df = counts_df.set_index(['organism', 'type'], append=True)
    counts_df = counts_df.reorder_levels(['organism', 'type', 'long_ID'])

    if feature_types_to_keep:
        counts_df = counts_df[counts_df.index.get_level_values('type').isin(feature_types_to_keep)]
    
    feature_df = counts_df.groupby(['type']).sum() 
    metadata_df = pd.concat([feature_df, metadata_df])

    return metadata_df, counts_df

def get_dge_table(results_table, attributes_df=None, fdr=0.1, sort=True, dge_out_path=None, rearrange_columns=True):
        if fdr:
            dge_table_subset = results_table[results_table['padj'] < fdr]
        else:
            dge_table_subset = results_table
        dge_table_subset = dge_table_subset.join(attributes_df)
        # Symmetric Log function (never undefined)
        def symlog10(x):
            if x > 0:
                return math.log10(x+1)
            elif x < 0:
                return -math.log10(-x+1)
            else:
                return 0
            
        dge_table_subset['symlog10baseMean'] = dge_table_subset['baseMean'].apply(lambda x: symlog10(x))

        if sort:
            dge_table_subset = dge_table_subset.sort_values('padj', ascending=True)
        if rearrange_columns:
            cols = list(dge_table_subset.columns.values)
            try:
                index_of_padj = cols.index('padj')
                index_of_baseMean = cols.index('baseMean')

                cols.remove('product')
                cols.insert(index_of_padj+1, 'product')

                dge_table_subset = dge_table_subset[cols]

                cols.remove('symlog10baseMean')
                cols.insert(index_of_baseMean+1, 'symlog10baseMean')

                dge_table_subset = dge_table_subset[cols]
            except ValueError:
                pass

        if dge_out_path:
            dge_table_subset.to_csv(dge_out_path, sep='\t')

        return dge_table_subset
    

def plot_feature_pct_sample(metadata_df, outdir, normalized=True):
    if normalized:
        read_counts = metadata_df.sum()
        metadata_df_pct_reads = metadata_df / read_counts
        axes = metadata_df_pct_reads.T.plot.bar(stacked=True)
        outpath = outdir / 'library_annotation_{}.png'.format('normalized')
    else:
        axes = metadata_df.T.plot.bar(stacked=True)
        outpath = outdir / 'library_annotation_{}.png'.format('sum')

    axes.legend(bbox_to_anchor=(1, 1))
    fig = axes.get_figure()
    fig.savefig(outpath, bbox_inches='tight', transparent=False, facecolor='w', edgecolor='w')


def plot_strain_featured_sample(counts_df, outdir, normalized = True):
    strain_counts_df = pd.DataFrame(columns=counts_df.columns)
    strain_counts_df = counts_df.copy()
    
    strain_counts_df = strain_counts_df.groupby('organism').sum()

    if normalized:
        strain_counts_df_pct_reads = strain_counts_df / strain_counts_df.sum()
        axes = strain_counts_df_pct_reads.T.plot.bar(stacked=True)
        outpath = outdir / 'strain_{}_featured_sample.png'.format('normalized')
    else:
        axes = strain_counts_df.T.plot.bar(stacked=True)
        outpath = outdir / 'strain_{}_featured_sample.png'.format('count')
    axes.legend(bbox_to_anchor=(1, 1))
    fig = axes.get_figure()
    fig.savefig(outpath, bbox_inches='tight', transparent=False, facecolor='w', edgecolor='w')

def plot_heatmap(count_table, indices=None, fig_out_path=None):
    log.info(f"plotting heatmap")
    log.debug(f"count_table:\n{count_table}")
    log.debug(f"indicies:\n{indices}")
    if indices is not None:
        log.debug(f"subsetted df:\n{count_table.loc[indices]}")
        figure = sns.clustermap(count_table.loc[indices])
    else:
        figure = sns.clustermap(count_table)
    if fig_out_path:
        figure.savefig(fig_out_path)

def get_PCA(count_table):
    pca = PCA(n_components=2)
    pca.fit(count_table)

    log.info(f'pca components: {pca.components_}')
    log.info(f'pca explained variance {pca.explained_variance_}')

    return pca

def plot_PCA(count_table, condition_table, fig_out_path):
    pca = PCA(n_components=2)
    log.info('DF Components')
    log.info(pd.DataFrame(pca.components_, columns = count_table.columns))
    log.info('PCA explained variance')
    explained_variance = list(pca.explained_variance_)
    log.info(explained_variance)
    transform = pca.fit_transform(count_table.T)
    
    samples = count_table.columns
    treatments = condition_table['treatment']
    unique_treatments = treatments.unique()

    color_dict = {}
    color=iter(plt.cm.rainbow(np.linspace(0,1,len(unique_treatments))))
    for i, t in enumerate(unique_treatments):
        color_dict[t] = next(color)

    fig = plt.figure(figsize=(10, 5), constrained_layout=True)
    ax = fig.gca()
    ax.scatter(transform[:,0], transform[:,1], c=list(map(color_dict.get, treatments)))

    for i, txt in enumerate(samples):
        ax.annotate(txt, transform[i])
    
    ax.set_xlabel('PC 1 (explained variance: {})'.format(explained_variance[0]))
    ax.set_ylabel('PC 2 (explained variance: {})'.format(explained_variance[1]))

    patches = [mpatches.Patch(color=c, label=l) for l, c in color_dict.items()]
    ax.legend(handles=patches)

    plt.savefig(fig_out_path)
    plt.close()

def plot_pct_genes_mapped_per_organism_sample(pct_genes_mapped_per_organism_sample_df, out_dir):
    for organism in pct_genes_mapped_per_organism_sample_df.index.unique():
        fig_out_path = out_dir / 'pct_genes_mapped_per_sample_{}.png'.format(organism)

        fig = plt.figure(figsize=(10, 5), constrained_layout=True)
        pct_genes_mapped_per_organism_sample_df.loc[organism].plot.bar()

        plt.savefig(fig_out_path)
        plt.close()

def plot_lfc_mean_normalized_counts(results_df, fig_out_path, fdr=0.1):
    
    fig = plt.figure(figsize=(10, 5), constrained_layout=True)
    ax = fig.gca()

    color_dict = {True : 'royalblue', False : 'darkgray'}
    color_map = list(map(color_dict.get, results_df['padj'] < fdr))

    ax.scatter(results_df['baseMean'], results_df['log2FoldChange'], c=color_map, alpha=0.5)
    
    ax.set_xscale('log')
    ax.set_xlabel('baseMean')
    ax.set_ylabel('log2FoldChange')
    
    plt.savefig(fig_out_path)
    plt.close()

def main(results_dir, htseq2dir, gff_dir, condition_table_path, raw_reads_dir, feature_types_to_keep=None):

    comparisons_dir = makeOutDir(results_dir, 'comparisons')
    metadata_figures_dir = makeOutDir(results_dir, 'metadata_figures')
    metadata_tables_dir = makeOutDir(results_dir, 'metadata_tables')

    app_dfs_dir = makeOutDir(results_dir, 'app_dfs')
    shutil.copy(condition_table_path, app_dfs_dir / condition_table_path.name)
    
    # attributes and annotations

    gffs = []
    for gff_path in gff_dir.iterdir():
        organism_name = str(gff_path.stem)
        annotation = gffpd.read_gff3(gff_path)
        attributes_df = annotation.attributes_to_columns()
        attributes_df['organism'] = organism_name
        gffs.append(attributes_df)

    attributes_df = pd.concat(gffs)

    attributes_df = attributes_df.rename(columns={'ID':'long_ID'})

    attributes_df = attributes_df.set_index('long_ID')

    organism_type_ID_df = attributes_df[['organism', 'type']]

    organism_type_ID_df.to_csv(results_dir / "organism_type_ID_df.tsv", sep="\t")

    metadata_df, counts_df = get_htseq2_and_metadata_df(htseq2dir, organism_type_ID_df, feature_types_to_keep)
    
    counts_df.to_csv(results_dir / "counts.tsv", sep="\t")
    counts_df.to_csv(app_dfs_dir / "counts.tsv", sep="\t")
    attributes_df.to_csv(results_dir / "attributes.tsv", sep="\t")

    metadata_df.to_csv(metadata_tables_dir / 'metadata.tsv', sep='\t')
    metadata_df.to_csv(app_dfs_dir / 'metadata.tsv', sep='\t')
    
    plot_feature_pct_sample(metadata_df, metadata_figures_dir, normalized=True)
    plot_feature_pct_sample(metadata_df, metadata_figures_dir, normalized=False)
    plot_strain_featured_sample(counts_df, metadata_figures_dir, normalized=True)
    plot_strain_featured_sample(counts_df, metadata_figures_dir, normalized=False)

    # pct of genes with >1 read for each organism
    organism_sample_gene_mapping_pct_df = pd.DataFrame(columns=counts_df.columns)
    count_df_gene_sampled_bool = counts_df >= 1
    count_df_gene_sampled_all = counts_df >= 0

    organism_sample_gene_mapping_pct_df = count_df_gene_sampled_bool.groupby('organism').sum() / count_df_gene_sampled_all.groupby('organism').sum()

    organism_sample_gene_mapping_pct_path = metadata_tables_dir / 'pct_genes_mapped_per_organism_sample.tsv'
    organism_sample_gene_mapping_pct_df.to_csv(organism_sample_gene_mapping_pct_path, sep='\t')

    plot_pct_genes_mapped_per_organism_sample(organism_sample_gene_mapping_pct_df, metadata_figures_dir)

    # comparisons and DEseq2
    conditions_df =  pd.read_csv(condition_table_path, sep='\t', index_col='sample_name')

    # extracts columns that include the word "comparison"
    comparisons_include_df = conditions_df[conditions_df.columns[np.array(conditions_df.columns.str.contains('comparison'))]]
    
    # creates new dataframe from comparisons_include_df that replaces control/treatment with True in order to subset future datasets
    conditions_df_sans_comparisons = conditions_df[conditions_df.columns[~np.array(conditions_df.columns.str.contains('comparison'))]]

    results_dfs = []
    rlog_dfs = []
    

    for comparison in comparisons_include_df.columns:
        
        log.info(comparison)
        
        comparison_dir = makeOutDir(comparisons_dir, comparison)
        comparison_count_dir = makeOutDir(comparison_dir, 'count_table')
        comparison_condition_dir = makeOutDir(comparison_dir, 'condition_table')
        comparison_deseq2out_dir = makeOutDir(comparison_dir, 'DEseq2out')
        figure_out_dir = makeOutDir(comparison_dir, 'post_DE_seq_figures')
        dge_out_dir = makeOutDir(comparison_dir, 'DGE_tables')

        included_samples_series = comparisons_include_df[comparisons_include_df[comparison].notna()][comparison]

        conditions = list(included_samples_series.unique())
        
        return_results_df = len(conditions) == 2 and 'control' in conditions and 'treatment' in conditions

        log.info(f"conditions {return_results_df} for {conditions}")

        comparison_condition_df = conditions_df_sans_comparisons.loc[included_samples_series.index]
        comparison_condition_df[comparison] = included_samples_series
        
        comparison_condition_df.to_csv(comparison_condition_dir / '{}_condition_table.tsv'.format(comparison), sep='\t')

        comparison_samples = included_samples_series.index.values 

        comparison_counts_df = counts_df[[f"sample_{s:02d}" for s in comparison_samples]]

        comparison_counts_df = comparison_counts_df.reset_index(level='type', drop=True)

        for organism in comparison_counts_df.index.unique(0):
            log.info(f'{organism}\n')

            comparison_organism_count_df = comparison_counts_df.loc[organism, :]

            # organism must be present in all samples. Good test to tell if it is from experience...
            # log.info('{} - organism lowest sample read mapping in comparison: {}'.format(organism, min(comparison_organism_count_df.mean())))
            if min(comparison_organism_count_df.mean()) > 1:
                
                comparison_raw_count_outpath = comparison_count_dir / f"{comparison}_raw_counts_{organism}.tsv"
                comparison_organism_count_df.to_csv(comparison_raw_count_outpath, sep='\t')

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

                    results_table_path = comparison_deseq2out_dir / '{}_{}_results.tsv'.format(comparison, organism)
                    results_table.to_csv(results_table_path, sep='\t')

                normalized_counts_df = pd.DataFrame(normalized_counts_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                rlog_df = pd.DataFrame(rlog_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                vst_df = pd.DataFrame(vst_array, index=comparison_organism_count_df.index, columns=included_samples_series.index.values)
                
                # saving out dataframes:
                rlog_count_data_path = comparison_deseq2out_dir / '{}_{}_rlog.tsv'.format(comparison, organism)
                vst_data_path = comparison_deseq2out_dir / '{}_{}_vst.tsv'.format(comparison, organism)
                normalized_counts_data_path = comparison_deseq2out_dir / '{}_{}_normalized_counts.tsv'.format(comparison, organism)
                
                normalized_counts_df.to_csv(normalized_counts_data_path, sep='\t')
                rlog_df.to_csv(rlog_count_data_path, sep='\t')
                vst_df.to_csv(vst_data_path, sep='\t')

                # post-DEseq2 analysis
                zero_dge_out_path = dge_out_dir / '{}_{}_DGE_all.tsv'.format(comparison, organism)
                dge_out_path = dge_out_dir / '{}_{}_DGE_fdr10pct.tsv'.format(comparison, organism)
                high_dge_out_path = dge_out_dir / '{}_{}_DGE_fdr01pct.tsv'.format(comparison, organism)
                pca_out_path = figure_out_dir / '{}_{}_PCA.png'.format(comparison, organism)
                clustermap_log_out_path = figure_out_dir / '{}_{}_clustermap_rlog.png'.format(comparison, organism)
                lfc_mean_normalized_counts_path = figure_out_dir / '{}_{}_lfc_v_mean_normalized_counts.png'.format(comparison, organism)

                plot_PCA(rlog_df, comparison_condition_df, pca_out_path)

                if return_results_df:

                    all_dge_table = get_dge_table(results_table, attributes_df, fdr=None, sort=True, dge_out_path=zero_dge_out_path)
                    mid_dge_table = get_dge_table(results_table, attributes_df, fdr=0.1, sort=True, dge_out_path=dge_out_path)

                    plot_lfc_mean_normalized_counts(results_table, lfc_mean_normalized_counts_path, fdr=0.1)

                    # get highest differentially expressed genes
                    high_dge_table = get_dge_table(results_table, attributes_df, fdr=0.01, sort=True, dge_out_path=high_dge_out_path)
                    try:
                        filtered_indicies = high_dge_table.index.tolist()
                        if len(filtered_indicies) == 0:
                            filtered_indicies = None
                        plot_heatmap(rlog_df, indices=filtered_indicies, fig_out_path=clustermap_log_out_path)
                    except:
                        log.error(f"rlog min: (something went wrong)\n{rlog_df.min()}")
                        raise
                    
                    #results df
                    index = pd.MultiIndex.from_product([[organism], results_table.index])
                    columns = pd.MultiIndex.from_product([[organism], [comparison], results_table.columns])
                    multiindex_results_df = results_df
                    multiindex_results_df.columns = columns
                    multiindex_results_df.index = index
                    results_dfs.append(multiindex_results_df)

                    #rlog df
                    index = pd.MultiIndex.from_product([[organism], rlog_df.index])
                    columns = pd.MultiIndex.from_product([[organism], [comparison], rlog_df.columns])
                    multiindex_rlog_df = rlog_df
                    multiindex_rlog_df.columns = columns
                    multiindex_rlog_df.index = index
                    rlog_dfs.append(multiindex_rlog_df)


