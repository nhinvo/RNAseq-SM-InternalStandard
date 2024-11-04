"""
This script performs filtering of internal standards. 

(more desc.)

Method Reference: S.M. Gifford (2016) "Quantitative Transcriptomics Reveals the Growth- and NutrientDependent Response 
of a Streamlined Marine Methylotroph to Methanol and Naturally Occurring Dissolved Organic Matter"
Link: https://journals.asm.org/doi/10.1128/mbio.01279-16

Molecular Weight calculation referece: "DNA and RNA Molecular Weights and Conversions", ThermoFisher
Link: https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html

Author: Nhi N. Vo 
Create Date: 11/15/23
"""

import pandas as pd 
import numpy as np
from Bio import SeqIO
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt

####################################################################################################
### INPUT PATHS ### 
sample_tsv_path = snakemake.config['input']['sample table']  # path to samples tsv file 
seq_path = snakemake.config['input']["reference genomes"]["InternalStandard"]["genome"]  # path to fasta of internal standard sequence
concentration_path = snakemake.config["input"]["internal_standard_concentration"]  # path to concentration added for each IS
uncorrected_counts_path = snakemake.input['uncorrected_counts']  # uncorrected counts of mapped reads 
annotation_path = snakemake.input['annotation_df']  # feature annotation (gff)

READ_LEN = snakemake.config["read length"]  # length of a sequenced read
STANDARD_SAMPLE_MIN = snakemake.config["minimum standard sample"]  # minimum number of sample a standard needs to be in 

### OUTPUT PATHS ###
standard_metadata_path = snakemake.output['standard_metadata']  # metadata tsv output path 
filtered_standard_metadata = snakemake.output['filtered_standard_metadata']  
plot_outdir = snakemake.output['IS_analysis_plot_outdir']  # path to dir for plots 
efficiency_slopes_intercepts = snakemake.output['efficiency_slopes_intercepts']  
####################################################################################################


def import_count_df():
    """
    Returns uncorrected count df with 3 columns: 
        - ['standard_name', 'sample', 'standard_recovered']
    """
    df = pd.read_csv(uncorrected_counts_path, sep='\t')
    df = df[df['organism'] == 'InternalStandard']  # counts for IS rows only 
    df = df.drop(['organism'], axis=1)

    # wide to long 
    df = df.melt(
        id_vars=['ID'], 
        var_name='sample', 
        value_name='reads_mapped'
    )

    df = df.rename(columns={'ID':'standard_name'})

    return df

def import_annot_df():
    """
    Returns Internal Standard length annotation df with 2 columns: 
        - ['standard_name', 'len']
    """
    df = pd.read_csv(annotation_path, sep='\t')
    df = df[df['organism'] == 'InternalStandard']  # counts for IS rows only 
    df['len'] = df['end'] - df['start']
    df = df[['ID', 'len']]
    df = df.rename(columns={'ID': 'standard_name'})

    return df 

def import_concentration_df():
    """
    Import df with internal standard concentration data and 
    check to make sure that all the required columns are 
    in place and in the correct format. 
    """
    df = pd.read_csv(concentration_path, sep='\t')
    required_columns = ['standard_name', 'standard_group', 'concentration (ng/ul)', 'volume_added (ul)']

    # Check if all required colunms are present in df
    if not df.columns.isin(required_columns).all():
        missing_columns = set(required_columns) - set(df.columns)
        error_message = f"Missing column(s) in concentration_df: {missing_columns}. Or column names aren't in the correct format."
        raise ValueError(error_message)

    # ensure dtypes are correct
    df = df.astype({
        'concentration (ng/ul)': float, 
        'volume_added (ul)': float, 
    })

    return df


def import_tables():
    """
    Returns imported pandas DataFrames from paths listed. 
    Ensures that the units are in correct format: 
        - for concentration_df: 
            - concentration column: 'concentration (ng/ul)'
            - volume column: 'volume_added (ul)'
    """
    # 1. sample tsv - from snakemake input directory
    sample_df = pd.read_csv(sample_tsv_path, sep='\t')

    # 2. uncorrected read mapping count - from mapping result 
    count_df = import_count_df()

    # 3. annotation df - for internal standard length 
    annot_df = import_annot_df()
    
    # 4. standard concentration - from snakemake input directory
    concentration_df = import_concentration_df()

    return sample_df, count_df, annot_df, concentration_df

def obtain_standard_recovered(count_df, annot_df):
    """
    Returns a DataFrame with number of standard RNA molecule (transcript) recovered
    with columns: ['standard_name', 'sample', 'standard_recovered']

    Number of transcript (molecule) of standard calculated by: 
        + (number of mapped reads * read len) / internal standard length
    """
    df = pd.merge(count_df, annot_df, on='standard_name', how='left')
    df['mapped_base'] = df['reads_mapped'] * READ_LEN
    df['standard_recovered'] = df['mapped_base'] / df['len']

    return df 

def cal_standard_mw():
    """
    Returns a pandas DataFrame of the molecular weight (mw) and length (seq_len) for 
    each internal standard added - calculated from the fasta sequence file. 
    """
    mol_data_list = []  

    for record in SeqIO.parse(seq_path, 'fasta'): 
        standard_name = record.id
        seq = record.seq.transcribe()  # obtain RNA from DNA seq
        length = len(seq)

        # calculate molecular weight based on sequence 
        mw_A = seq.count('A') * 329.2
        mw_U = seq.count('U') * 306.2
        mw_C = seq.count('C') * 305.2
        mw_G = seq.count('G') * 345.2
        constant = 159  # note: addition of "159" = takes into account the M.W. of a 5' triphosphate.
        mw = mw_A + mw_U + mw_C + mw_C + constant 

        mol_data_dict = {
            'standard_name': standard_name,  # name should match with concentration_df
            'mw': float(mw), 
            'seq_len': int(length)
        }

        mol_data_list.append(mol_data_dict)

    mw_df = pd.DataFrame(mol_data_list)

    return mw_df

def obtain_standard_added(concentration_df, mw_df):
    """
    Returns a DataFrame with number of standard molecule added - calculated from 
    standard molecular weight, mass, and volume. 
    """
    # combine mw with concentration & volume 
    added_df = concentration_df.join(mw_df.set_index('standard_name'), on='standard_name')

    # calculate mass from concentration & volume 
    added_df['mass_ng'] = added_df['concentration (ng/ul)'] * added_df['volume_added (ul)'] 
    added_df['mass_g'] = added_df['mass_ng'] / 1000000000  # ng to g

    # calculate number of molecule (standard added) from mass and mw
    added_df['moles_added'] = added_df['mass_g'] / added_df['mw']
    added_df['standard_added'] = (
        (added_df['moles_added'] * (6.022 * (10**23)))
        .astype(int)  # round to whole number of standard molecule added 
    )

    return added_df 

def obtain_standard_data(sample_df, count_df, annot_df, concentration_df):
    """
    Obtains internal standard metadata, such as: 
        - Number of molecules (transcript) recovered ((reads mapped * read len) / internal standard length) 
        - Number of molecules (transcript) added based on: concentration, volume, seq length
    """    
    # calculate number of standard transcript added 
    recovered_df = obtain_standard_recovered(count_df, annot_df)

    # calculate molecular weight of standard RNA molecule added 
    mw_df = cal_standard_mw()
    added_df = obtain_standard_added(concentration_df, mw_df)
    
    # combine added and recovered data together 
    standard_df = added_df.join(recovered_df.set_index('standard_name'), on=['standard_name'])

    # logging
    standard_df = standard_df[standard_df['standard_recovered']!=0]  # remove standards where no reads were recovered 
    standard_df['standard_added_log10'] = np.log10(standard_df['standard_added'])
    standard_df['standard_recovered_log10'] = np.log10(standard_df['standard_recovered'])

    # filter to remove standards recovered in less than specified number of samples
    standard_groups = standard_df.groupby(['standard_name'])
    filtered_dfs = []
    for index, gdf in standard_groups:
        if len(gdf) < STANDARD_SAMPLE_MIN:
            # skip standards with less than specified number of samples 
            print(f"Removing standard {index[0]}; reason: recovered in less than {STANDARD_SAMPLE_MIN} samples.")
            continue
        else: 
            filtered_dfs.append(gdf)
    # recombine standards that pass sample filter 
    standard_df = pd.concat(filtered_dfs)

    # add color for downstream plotting 
    iwanthue = [
        "#da8245",
        "#9f92ff",
        "#58cc9b",
        "#5c3c90",
        "#5a9a41",
        "#be59a2",
        "#c1ba52",
        "#ca4e69",
        "#8d7524",
        "#a33a2b"
    ]

    colors = dict(zip(
        standard_df['standard_name'].unique(),
        iwanthue, 
    ))
    standard_df['color'] = standard_df['standard_name'].map(colors)

    standard_df.to_csv(standard_metadata_path, sep='\t', index=False)

    return standard_df

def obtain_linregress_data(df, y_col):
    """
    Performs linear regression on the dataframe provided. Calculates standard deviation, upper, lower bounds. 
    Returns regression results. 
    """
    # perform linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['standard_added_log10'], df[y_col])

    # standard deviation calculation 
    std_deviation = std_err * 3
    
    # Calculate upper and lower bounds by adding and subtracting the standard deviation from the regression line 
    upper_bound = (slope * df['standard_added_log10']) + intercept + std_deviation
    lower_bound = (slope * df['standard_added_log10']) + intercept - std_deviation

    x_range = np.linspace(5, 9.5, 100)  # extend x-values beyond current dataset 
    extended_upper_bound = (slope * x_range) + intercept + std_deviation
    extended_lower_bound = (slope * x_range) + intercept - std_deviation

    data={
        'df':df, 'slope':slope, 'intercept':intercept, 'upper_bound':upper_bound, 'lower_bound':lower_bound,
        'extended_upper_bound': extended_upper_bound, 'extended_lower_bound': extended_lower_bound, 
    }

    return data

def obtain_standard_median(df):
    """
    Returns df with columns:
        - standard_name: name of standard (same as that in input df)
        - median_standard_recovered_log10: for each standard, the median of 
            all standard_recovered_log10 values
        - standard_added_log10: value of standard added (same as that in input df)
        - color: standard color in plot (same as input df)
    """
    standard_groups = df.groupby(['standard_name'])

    data = {
        'standard_name': [], 
        'median_standard_recovered_log10': [], 
        'standard_added_log10': [], 
        'color': [], 
    }

    for index, gdf in standard_groups:
        data['standard_name'].append(index[0])
        data['median_standard_recovered_log10'].append(gdf['standard_recovered_log10'].median() )
        data['standard_added_log10'].append(gdf['standard_added_log10'].iloc[0])
        data['color'].append(gdf['color'].iloc[0])

    df = pd.DataFrame(data)

    return df

def plot_all_iterations(data_dict_list):
    """
    For each median filtering iteration, plot median internal standards.
    """
    # make sure plot outdir is there 
    linregress_plot_outdir = Path(f"{plot_outdir}/linregress_plots_median")
    linregress_plot_outdir.mkdir(parents=True, exist_ok=True) 

    for iteration, data_dict in enumerate(data_dict_list):
        df = data_dict['df']

        fig, ax = plt.subplots(1, 1, figsize=(10,7))

        # scatter plotting 
        plt.scatter(df['standard_added_log10'], df['median_standard_recovered_log10'], s=120, zorder=2, facecolors='none', edgecolors=df['color'], linewidths=3)

        # regression line + upper & lower bounds 
        x_range = np.linspace(5, 9.5, 100)  # extend x-values beyond current dataset 
        plt.plot(x_range, (data_dict['slope'] * x_range) + data_dict['intercept'], color='tab:purple', zorder=1, linewidth=2)
        plt.plot(x_range, data_dict['extended_upper_bound'], color='tab:gray', label='Upper Bound', zorder=1, linewidth=2)
        plt.plot(x_range, data_dict['extended_lower_bound'], color='tab:gray', label='Lower Bound', zorder=1, linewidth=2)

        # axis limits
        ax.set_ylim([-2, 5])
        ax.set_xlim([5, 9.5])

        fig.supxlabel('log10(standard_added)')
        fig.supylabel('median log10(reads_recovered)')

        plt.tight_layout()
        plt.savefig(f"{linregress_plot_outdir}/median_standard_added_vs_recovered_{iteration+1}.png")

        plt.close()  

def IS_filtering_median(df):
    """
    Returns a list of dictionary storing data for before and after filtering. 
    """
    # obtain median for each standard 
    df = obtain_standard_median(df)

    # plotting data
    all_linregress_data = [] 

    # obtain lingress data before filtering
    linregress_data = obtain_linregress_data(df, 'median_standard_recovered_log10')
    all_linregress_data.append(linregress_data)  # save that data for plotting 

    # test to see if any outliers outside bounds 
    upper_bound_condition = df['median_standard_recovered_log10'] > linregress_data['upper_bound']
    lower_bound_condition = df['median_standard_recovered_log10'] < linregress_data['lower_bound']
    outliers = upper_bound_condition | lower_bound_condition

    if outliers.any():  # if there were outliers based on conditions above 
        df = df[~outliers]  # remove outliers from sample_df

        # obtain lingress data 
        linregress_data = obtain_linregress_data(df, 'median_standard_recovered_log10')
        all_linregress_data.append(linregress_data)  # save that data for plotting 

    # plot data for each global filtering iteration 
    plot_all_iterations(all_linregress_data)

    return all_linregress_data

def plot_all_samples(df, slope_int_df, fname_count):
    """
    For each sample in df, plot internal standards. 
    """
    # make sure plot outdir is there 
    sample_plot_outdir = Path(f"{plot_outdir}/linregress_plots_all_samples")
    sample_plot_outdir.mkdir(parents=True, exist_ok=True) 

    # group by samples
    sample_groups = df.groupby(['sample'])
    sample_groups = [(index[0], sdf) for index, sdf in sample_groups]

    # values for figure setup 
    num_subplots = len(sample_groups)
    num_cols = int(num_subplots ** 0.5 + 1)
    num_rows = (num_subplots + num_cols - 1) // num_cols

    # make figure
    figsize = (num_cols * 4.5, num_rows * 4.5)
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)

    for count, ax in enumerate(fig.axes):
        try:
            sample_data = sample_groups[count]
        except IndexError:
            "More subplots specified than samples being plotted. Leaving subplot empty."
            ax.set_axis_off()
            continue

        sample = sample_data[0]
        sdf = sample_data[1]
        
        slope = slope_int_df[slope_int_df['sample'] == sample]['slope'].iloc[0]
        intercept = slope_int_df[slope_int_df['sample'] == sample]['intercept'].iloc[0]

        # scatter plotting 
        ax.scatter(sdf['standard_added_log10'], sdf['standard_recovered_log10'], s=120, zorder=2, facecolors='none', edgecolors=sdf['color'], linewidths=3)

        # regression plotting 
        x_range = np.linspace(5, 9.5, 100)  # extend x-values beyond current dataset 
        ax.plot(x_range, (slope * x_range) + intercept, color='tab:grey', zorder=1, linewidth=2)
        
        # axis limits
        ax.set_ylim([-2, 5])
        ax.set_xlim([5, 9.5])

        # titles
        ax.set_title(f"{sample} - slope={slope:.3f}")
        ax.set_xlabel('log10(standard_added)')
        ax.set_ylabel('log10(reads_recovered)')

    plt.subplots_adjust(top=0.995)  # makes it so that there is a little extra space on top of the plot - so subplot title doesn't get cutoff 
    plt.savefig(f"{sample_plot_outdir}/all_sample_plots_{fname_count+1}.png")

    plt.close()  

def obtain_sample_slope_intercept(standard_dict_list, df):
    """
    For each filtering iteration (including before filtering), plot
    all internal standards and the linear fit. 
    """
    # plotting - every iteration 
    for count, iteration in enumerate(standard_dict_list):
        iter_standard_df = iteration['df']

        # list of standard names in df 
        standard_names = iter_standard_df['standard_name'].values.tolist()

        # filter standard df for standard in list 
        iter_df = df[df['standard_name'].isin(standard_names)].copy()

        # for each sample, obtain slope and intercept based on standards
        sample_groups = iter_df.groupby(['sample'])
        slope_int_data = {
            "sample": [], 
            "slope": [], 
            "intercept": [], 
        }
        for index, sdf in sample_groups:
            sname = index[0]
            lingress_data = obtain_linregress_data(sdf, 'standard_recovered_log10')

            slope_int_data['sample'].append(sname)
            slope_int_data['slope'].append(lingress_data['slope'])
            slope_int_data['intercept'].append(lingress_data['intercept'])

        
        slope_int_df = pd.DataFrame(slope_int_data)
        plot_all_samples(iter_df, slope_int_df, count)

        if count == (len(standard_dict_list) - 1):
            # last iteration (filtered data)
            slope_int_df.to_csv(efficiency_slopes_intercepts, sep='\t', index=False)
            iter_df.to_csv(filtered_standard_metadata, sep='\t', index=False)

    return slope_int_df

def plot_slope_int_distr(df):
    """
    """
    # make sure plot outdir is there 
    sample_plot_outdir = Path(f"{plot_outdir}/slope_int_distribution")
    sample_plot_outdir.mkdir(parents=True, exist_ok=True) 

    num_bins = int(np.sqrt(len(df)))
    # num_bins = int(np.log2(len(df)) + 1) 
    # num_bins = int(2 * len(df) ** (1/3))
    num_bins = 10

    cols = ['slope', 'intercept']
    for col in cols:
        fig, ax = plt.subplots()
        ax.hist(df[col], bins=num_bins)

        ax.set_title(f"Distribution of {col.title()}; bins={num_bins}")
        ax.set_ylabel('Sample Count')
        ax.set_xlabel(f'Internal Standard {col.title()}')

        plt.savefig(f"{sample_plot_outdir}/{col}.png")

        plt.close()  

if __name__ == "__main__":
    # import and process input tables 
    sample_df, count_df, annot_df, concentration_df = import_tables()

    # obtain standard data (molecule added vs. molecule recovered)
    standard_df = obtain_standard_data(sample_df, count_df, annot_df, concentration_df)

    # global filtering (obtain median of standard across all samples & filter)
    all_linregress_data = IS_filtering_median(standard_df)

    # for the remaining standards, obtain slope of those standards within each sample 
    final_slope_int_df = obtain_sample_slope_intercept(all_linregress_data, standard_df)

    # plot distribution of slope & intercept of filtered internal standards
    plot_slope_int_distr(final_slope_int_df)