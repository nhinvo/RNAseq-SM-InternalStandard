rule make_results_dir:
    output:
        directory(results_dir)
    shell:
        "sleep 1"

# run post-HTseq script
rule post_htseq2_parsing:
    input:
        sample_counts=expand(feature_count_dir / "{sample}.tsv", sample=SAMPLES),
        raw_gff_dir=gff_refs_dir,
        condition_table_path=condition_table_file,
        r_dir = results_dir,
        raw_reads=raw_reads_dir,
        raw_reads_counts=library_count_dir / "library_len.tsv"
    output:
        done_flag = touch(done_file_dir / "post_htseq2_parsing.done"),
    conda:
        "../envs/post_htseq2_parsing.yaml"
    script:
        "../scripts/post_htseq2_parsing_snakemake.py"
