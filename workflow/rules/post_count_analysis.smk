
# # run post-HTseq script
# rule generate_counts_metadata_dfs:
#     input:
#         sample_counts=expand(output_path_dict["feature_count"] / "{sample}.tsv", sample=SAMPLES),
#         raw_gff_dir=Path(config["input"]["gff_refs"]),
#         condition_table_path=Path(config["samples"]),
#         config_yaml_path=Path(config["config"]),
#         ref_json = results_path_dict['ref_json']
#     output:
#         counts = results_path_dict['counts'],
#         metadata = results_path_dict['metadata'],
#         config = results_path_dict['config'],
#         samples = results_path_dict['samples'],
#         counts_json = results_path_dict['counts_json'],
#         annotation_json = results_path_dict['annotation_json'],
#     conda:
#         "../envs/post_count_analysis.yaml"
#     benchmark:
#         benchmark_dir / 'gff_tools' / 'annotation_concat.benchmark'
#     resources:
#         partition = 'sched_mit_chisholm',
#         mem = '250G',
#         ntasks = 20,
#         time = '0-12', 
#         output = str(log_dir / 'gff_tools' / 'annotation_concat.out'),
#         error = str(log_dir / 'gff_tools' / 'annotation_concat.err'),
#     script:
#         "../scripts/generate_counts_metadata_dfs.py"

# # run post-HTseq script
# rule generate_ref_table:
#     input:
#         raw_gff_dir = Path(config["input"]["gff_refs"]),
#     output:
#         ref_json = results_path_dict['ref_json']
#     conda:
#         "../envs/post_count_analysis.yaml"
#     benchmark:
#         benchmark_dir / 'gff_tools' / 'annotation_concat.benchmark'
#     resources:
#         partition = 'sched_mit_chisholm',
#         mem = '250G',
#         ntasks = 20,
#         time = '0-12', 
#         output = str(log_dir / 'gff_tools' / 'annotation_concat.out'),
#         error = str(log_dir / 'gff_tools' / 'annotation_concat.err'),
#     script:
#         "../scripts/generate_bio_db_ref_table.py"

# rule generate_comparison_results:
#     input: 
#         comparison_yaml = Path(config["comparisons"]),
#         counts = results_path_dict['counts'],
#     output:
#         comparison_data = results_path_dict['comparisons']
#     conda:
#         "../envs/DEseq2py.yaml"
#     benchmark:
#         benchmark_dir / 'gff_tools' / 'annotation_concat.benchmark'
#     resources:
#         partition = 'sched_mit_chisholm',
#         mem = '250G',
#         ntasks = 20,
#         time = '0-12', 
#         output = str(log_dir / 'gff_tools' / 'annotation_concat.out'),
#         error = str(log_dir / 'gff_tools' / 'annotation_concat.err'),
#     script:
#         "../scripts/generate_comparison_results.py"

# rule generate_data_json:
#     input:
#         counts_json = results_path_dict['counts_json'],
#         annotation_json = results_path_dict['annotation_json'],
#         metadata = results_path_dict['metadata'],
#         config = results_path_dict['config'],
#         samples = results_path_dict['samples'],
#         comparison_data = results_path_dict['comparisons'],
#     output:
#         data_json = results_path_dict['json'],
#         samples = results_path_dict['samples_json'],
#         metadata = results_path_dict['metadata_json'],
#         config = results_path_dict['config_json'],
#     conda:
#         "../envs/post_count_analysis.yaml"
#     benchmark:
#         benchmark_dir / 'gff_tools' / 'annotation_concat.benchmark'
#     resources:
#         partition = 'sched_mit_chisholm',
#         mem = '250G',
#         ntasks = 20,
#         time = '0-12', 
#         output = str(log_dir / 'gff_tools' / 'annotation_concat.out'),
#         error = str(log_dir / 'gff_tools' / 'annotation_concat.err'),
#     script:
#         "../scripts/generate_data_json.py"