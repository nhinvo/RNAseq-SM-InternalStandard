
# run post-HTseq script
rule generate_counts_metadata_dfs:
    input:
        sample_counts=expand(output_path_dict["feature_count"] / "{sample}.tsv", sample=SAMPLES),
        raw_gff_dir=Path(config["input"]["gff_refs"]),
        condition_table_path=Path(config["samples"]),
        config_yaml_path=Path(config["config"]),
    output:
        counts = results_path_dict['counts'],
        metadata = results_path_dict['metadata'],
        config = results_path_dict['config'],
        samples = results_path_dict['samples'],
    conda:
        "../envs/post_htseq2_parsing.yaml"
    script:
        "../scripts/generate_counts_metadata_dfs.py"

# run post-HTseq script
rule generate_ref_table:
    input:
        raw_gff_dir = Path(config["input"]["gff_refs"]),
    output:
        ref_table = results_path_dict['ref_table']
    conda:
        "../envs/post_htseq2_parsing.yaml"
    script:
        "../scripts/generate_bio_db_ref_table.py"

rule generate_data_json:
    input:
        counts = results_path_dict['counts'],
        metadata = results_path_dict['metadata'],
        config = results_path_dict['config'],
        samples = results_path_dict['samples'],
        ref_table = results_path_dict['ref_table']
    output:
        data_json = results_path_dict['json']
    conda:
        "../envs/post_htseq2_parsing.yaml"
    script:
        "../scripts/generate_data_json.py"