rule make_results_dir:
    output:
        directory(config["results"])
    shell:
        "mkdir -p {output}; sleep 1"

# run post-HTseq script
rule post_htseq2_parsing:
    input:
        sample_counts=expand(path_dict["feature_count"] / "{sample}.tsv", sample=SAMPLES),
        raw_gff_dir=Path(config["input"]["gff_refs"]),
        condition_table_path=Path(config["samples"]),
        r_dir = Path(config["results"]),
        raw_reads=Path(config["input"]["raw_reads"]),
        raw_reads_counts=path_dict["library_count"],
        bam_coverages=expand(path_dict["coverage_positions"] / "{sample}_coverage_depth.tsv", sample=SAMPLES),
    output:
        done_flag = touch(path_dict["done_files"] / "post_htseq2_parsing.done"),
    conda:
        "../envs/post_htseq2_parsing.yaml"
    script:
        "../scripts/post_htseq2_parsing_snakemake.py"
