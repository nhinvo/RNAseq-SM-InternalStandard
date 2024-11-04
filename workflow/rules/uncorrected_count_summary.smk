rule generate_ref_table:
    output:
        results_dict['bio_db_ref']
    conda:
        "../envs/uncorrected_count_summary.yaml"
    script:
        "../scripts/generate_bio_db_ref_table.py"


rule generate_annotation_df:
    output:
        results_dict['annotations']
    conda:
        "../envs/uncorrected_count_summary.yaml"
    script:
        "../scripts/generate_annotation_df.py"


rule generate_raw_counts_metadata_dfs:
    """
    Obtain uncorrected counts (from htseq). 
    """
    input:
        annotations = results_dict['annotations'],
        counts = expand(scratch_dict["feature_count"] / "{sample}.tsv", sample=SAMPLES),
    output:
        raw_counts = results_dict['uncorrected_counts'],
        raw_counts_w_annotations = results_dict['raw_counts_w_annotations'],
        mapping_metadata = results_dict['mapping_metadata'],
        organism_occurance = results_dict['organism_occurance'],
        gene_sparsity = results_dict['gene_sparsity'],
    conda:
        "../envs/uncorrected_count_summary.yaml"
    script:
        "../scripts/generate_raw_counts_metadata_dfs.py"
