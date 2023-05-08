rule generate_ref_table:
    output:
        results_dict['bio_db_ref']
    conda:
        "../envs/post_count_analysis.yaml"
    script:
        "../scripts/generate_bio_db_ref_table.py"


rule generate_annotation_df:
    output:
        results_dict['annotations']
    conda:
        "../envs/post_count_analysis.yaml"
    script:
        "../scripts/generate_annotation_df.py"


rule generate_raw_counts_metadata_dfs:
    input:
        annotations = results_dict['annotations'],
        counts = expand(scratch_dict["feature_count"] / "{sample}.tsv", sample=SAMPLES),
    output:
        counts = results_dict['raw_counts'],
        counts_w_annotations = results_dict['counts_w_annotations'],
        mapping_metadata = results_dict['mapping_metadata'],
        organism_occurance = results_dict['organism_occurance'],
        gene_sparsity = results_dict['gene_sparsity'],
    conda:
        "../envs/post_count_analysis.yaml"
    script:
        "../scripts/generate_raw_counts_metadata_dfs.py"

rule generate_comparison_results:
    input: 
        counts = results_dict['raw_counts'],
        annotations = results_dict['annotations'],
    output:
        rlog = results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "rlog.tsv",
        vst = results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "vst.tsv",
        results = results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "results.tsv",
        results_w_annot = results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "results_w_annotation.tsv",
    conda:
        "../envs/DEseq2.yaml"
    script:
        "../scripts/generate_comparison_results.R"

rule generate_data_json:
    input:
        rlogs = expand(results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "rlog.tsv", zip, organism=ORGANISMS, experiment=EXPERIMENTS, comparison=COMPARISONS),
        vsts = expand(results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "vst.tsv", zip, organism=ORGANISMS, experiment=EXPERIMENTS, comparison=COMPARISONS),
        deseq2_results = expand(results_dict['DEseq2'] / "{organism}" / "{experiment}" / "{comparison}" / "results.tsv", zip, organism=ORGANISMS, experiment=EXPERIMENTS, comparison=COMPARISONS),
        annotations = results_dict['annotations'],
        counts = results_dict['raw_counts'],
        mapping_metadata = results_dict['mapping_metadata'],
        organism_occurance = results_dict['organism_occurance'],
        gene_sparsity = results_dict['gene_sparsity'],
        bio_df_ref = results_dict['bio_db_ref'],
    output:
        results_dict['data_json']
    conda:
        "../envs/post_count_analysis.yaml"
    script:
        "../scripts/generate_data_json.py"
