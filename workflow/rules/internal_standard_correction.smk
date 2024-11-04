rule internal_standard_analysis:
    """
    Performs internal standard (IS) filtering and obtain read sequencing efficiency. 
    For each sample: 
        - Obtain linear regression of all IS in each sample
        - Perform filtering to remove IS that don't fall within regression 
        - Save metadata of each IS
        - Obtain slope and intercept of linear regression from IS 
            - For back-calculating (correction) of gene read count 
    """
    input:
        uncorrected_counts = results_dict['uncorrected_counts'],
        annotation_df = results_dict['annotations']
    output:
        # metadata (len, mw, num_added, etc.) for the all standards added 
        standard_metadata=results_dict['all_standard_metadata'],  

        # metadata of standard after filtering out outliers 
        filtered_standard_metadata=results_dict['filtered_standard_metadata'],  

        # directory for IS-filtering plots 
        IS_analysis_plot_outdir = directory(results_dict['IS_analysis_plot_outdir']),
         
        # slope and intercept from standards regression of each sample 
        efficiency_slopes_intercepts=results_dict['efficiency_slopes_intercepts'], 
    conda:
        "../envs/internal_standard_analysis.yaml"
    script:
        "../scripts/internal_standard_analysis.py"


rule transcript_back_calculate:
    """
    Performs read count correction on raw mapped read counts (uncorrected count) using
    slopes and intercept obtained from internal standard filtering/analysis. 
    """
    input:
        uncorrected_counts = results_dict['uncorrected_counts'], 
        efficiency_slopes_intercepts = results_dict['efficiency_slopes_intercepts'], 
        annotation_df = results_dict['annotations'],
    output:
        absolute_counts = results_dict['absolute_gene_count_per_cell'],
        summary_back_calculate_counts = results_dict['summary_back_calculate_counts'],
        corrected_transcript_count = results_dict['corrected_transcript_count'], 
    conda:
        "../envs/internal_standard_analysis.yaml"
    script:
        "../scripts/transcript_back_calculate.py"
