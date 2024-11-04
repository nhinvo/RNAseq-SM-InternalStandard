""" 
This script calculates transcripts/library and transcript/cell using slope and
intercept of the linear line obtained from filtered internal standards. Then 
adjust the corrected transcript count by gene length. 

Author: Nhi N. Vo 
Date: 12/03/23
"""

import math
import pandas as pd 
import numpy as np
from pathlib import Path

####################################################################################################
### INPUT PATHS ###
standard_slope_int = snakemake.input['efficiency_slopes_intercepts']
uncorrected_counts_path = snakemake.input['uncorrected_counts']  # uncorrected count data (raw reads mapped)
cell_count_path = snakemake.config["input"]["cell_count_data"]  # path to cell count data 
annotation_path = snakemake.input['annotation_df']  # feature annotation (gff)

READ_LEN = snakemake.config["read length"]  # length of a sequenced read

### OUTPUT PATHS ###
corrected_count_outpath = snakemake.output['corrected_transcript_count'] 
summary_back_calculate_counts_outpath = snakemake.output['summary_back_calculate_counts'] 
absolute_gene_count_per_cell_outpath = snakemake.output['absolute_counts']
####################################################################################################


def import_annot_df():
    """
    Returns length annotation (for all genes) df with 2 columns: 
        - ['ID', 'gene_len']
    """
    df = pd.read_table(annotation_path)
    df = df[df['organism']!='InternalStandard']  # remove IS count rows 
    df['gene_len'] = df['end'] - df['start']
    df = df[['ID', 'gene_len']]

    return df

def process_read_count_df(annot_df):
    """
    Returns df with columns: ['organism', 'ID', 'sample', 'mapped_reads', 'gene_len', 'transcript_count']
        - Obtain number of transcript molecule recovered (mapped) by: 
            + (number of mapped reads * read len) / internal standard length
    """
    df = pd.read_table(uncorrected_counts_path)
    df = df[df['organism']!='InternalStandard']  # remove IS count rows  

    # transform table so cols are: [organism, ID, sample, count]
    df = df.melt(id_vars=['organism', 'ID'], var_name='sample', value_name='mapped_reads')

    df = pd.merge(df, annot_df, on='ID', how='left')

    df['mapped_base'] = df['mapped_reads'] * READ_LEN
    df['transcript_count'] = df['mapped_base'] / df['gene_len']

    df = df.drop(columns=['mapped_base'])

    return df

def process_cell_count_df():
    """
    """
    cell_count_df = pd.read_csv(cell_count_path, sep='\t')

    try:
        cell_count_df = cell_count_df[['sample', 'total_cell_count']]
    except:
        error_message = f"Missing column(s) in cell_count_df, or column names aren't in the correct format. Required columns are: ['sample', 'total_cell_count']"
        raise ValueError(error_message)

    no_cell_count_samples = cell_count_df[cell_count_df['total_cell_count'].isnull()]['sample'].values.tolist()
    print(f"Samples without cell count data: {no_cell_count_samples}.")

    # remove samples without count data
    cell_count_df = cell_count_df[cell_count_df['total_cell_count'].notna()]  

    # ensure dtypes are correct
    cell_count_df = cell_count_df.astype({
        'total_cell_count': int, 
    })

    return cell_count_df

def import_dfs():
    """
    Import and process dataframes for downstream back-calculation. 
    """
    # import standard slope and intercept data 
    standard_slope_int_df = pd.read_table(standard_slope_int, sep='\t')

    # import annotation df (gene length data)
    annot_df = import_annot_df()

    # import read count data and calculate number of transcript molecule recovered
    transcript_count_df = process_read_count_df(annot_df)

    # import cell count data 
    cell_count_df = process_cell_count_df()

    return standard_slope_int_df, transcript_count_df, cell_count_df

def correct_transcript_count(y, m, b):
    """
    y=reads sequenced; m=slope; b=intercept
    """
    try:
        x = 10**((math.log10(y)-b)/m)
    except:  # for when y (reads sequenced) is 0 
        x = False
    return x

def back_calculate(standard_slope_int_df, transcript_count_df, cell_count_df):
    """
    """
    sdfs = []
    sample_groups = transcript_count_df.groupby(['sample'])

    for sample_name, sdf in sample_groups:
        sname = sample_name[0]
        
        # obtain slope and intercept of internal standard linear line 
        m = standard_slope_int_df[standard_slope_int_df['sample']==sname]['slope'].iloc[0]
        b = standard_slope_int_df[standard_slope_int_df['sample']==sname]['intercept'].iloc[0]

        # obtain corrected count using transcript added and IS effeciency (x from y): x=10^((log10(y)-b)/m)
        sdf['corrected_transcript_count'] = sdf.apply(lambda row : correct_transcript_count(row['transcript_count'], m, b), axis = 1)

        # obtain transcript/cell - make rows empty for samples without cell count data 
        if sname not in cell_count_df['sample'].values.tolist():
            sdf['count_per_cell'] = np.nan  
            sdfs.append(sdf)
            continue 

        sample_cell_count = cell_count_df[cell_count_df['sample']==sname]['total_cell_count'].iloc[0]
        sdf['count_per_cell'] = sdf['corrected_transcript_count'].div(sample_cell_count)
        sdfs.append(sdf)

    count_per_cell_df = pd.concat(sdfs)
    corrected_count_df = count_per_cell_df[['organism', 'sample', 'ID', 'corrected_transcript_count']].copy()
    gene_count_per_cell_df = count_per_cell_df[['organism', 'sample', 'ID', 'count_per_cell']].copy()

    corrected_count_df = corrected_count_df.pivot(
        index=['organism', 'ID'], columns='sample', values='corrected_transcript_count'
    ).reset_index()

    gene_count_per_cell_df = gene_count_per_cell_df.pivot(
        index=['organism', 'ID'], columns='sample', values='count_per_cell'
    ).reset_index()

    count_per_cell_df.to_csv(f"{summary_back_calculate_counts_outpath}", sep='\t', index=False) 
    corrected_count_df.to_csv(f"{corrected_count_outpath}", sep='\t', index=False) 
    gene_count_per_cell_df.to_csv(f"{absolute_gene_count_per_cell_outpath}", sep='\t', index=False) 


if __name__ == "__main__":
    standard_slope_int_df, transcript_count_df, cell_count_df = import_dfs()

    back_calculate(standard_slope_int_df, transcript_count_df, cell_count_df)