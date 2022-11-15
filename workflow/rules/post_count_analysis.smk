rule generate_ref_table:
    output:
        results_dict['bio_db_ref']
    conda:
        "../envs/post_count_analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '1G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'post_count_analysis' / 'generate_ref_table.out'),
        error = lambda wildcards: mk_err(log_dir / 'post_count_analysis' / 'generate_ref_table.err'),
    script:
        "../scripts/generate_bio_db_ref_table.py"


rule generate_annotation_df:
    output:
        results_dict['annotations']
    conda:
        "../envs/post_count_analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'post_count_analysis' / 'generate_annotation_df.out'),
        error = lambda wildcards: mk_err(log_dir / 'post_count_analysis' / 'generate_annotation_df.err'),
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
    resources:
        partition = 'sched_mit_chisholm',
        mem = '1G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'post_count_analysis' / 'generate_raw_counts_metadata_dfs.out'),
        error = lambda wildcards: mk_err(log_dir / 'post_count_analysis' / 'generate_raw_counts_metadata_dfs.err'),
    script:
        "../scripts/generate_raw_counts_metadata_dfs.py"

rule generate_comparison_results:
    input: 
        counts = results_dict['raw_counts'],
        annotations = results_dict['annotations'],
    output:
        rlog = results_dict['DEseq2'] / "{organism}" / "{comparison}" / "rlog.tsv",
        vst = results_dict['DEseq2'] / "{organism}" / "{comparison}" / "vst.tsv",
        results = results_dict['DEseq2'] / "{organism}" / "{comparison}" / "results.tsv",
        results_w_annot = results_dict['DEseq2'] / "{organism}" / "{comparison}" / "results_w_annotation.tsv",
    conda:
        "../envs/DEseq2py.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '1G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'post_count_analysis', f'{wildcards.organism} {wildcards.comparison}'),
        error = lambda wildcards: mk_err(log_dir / 'post_count_analysis', f'{wildcards.organism} {wildcards.comparison}'),
    script:
        "../scripts/generate_comparison_results.py"

rule generate_data_json:
    input:
        rlogs = expand(results_dict['DEseq2'] / "{organism}" / "{comparison}" / "rlog.tsv", zip, organism=ORGANISMS, comparison=COMPARISONS),
        vsts = expand(results_dict['DEseq2'] / "{organism}" / "{comparison}" / "vst.tsv", zip, organism=ORGANISMS, comparison=COMPARISONS),
        deseq2_results = expand(results_dict['DEseq2'] / "{organism}" / "{comparison}" / "results.tsv", zip, organism=ORGANISMS, comparison=COMPARISONS),
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
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'post_count_analysis' / 'generate_data_json.out'),
        error = lambda wildcards: mk_err(log_dir / 'post_count_analysis' / 'generate_data_json.err'),
    script:
        "../scripts/generate_data_json.py"
