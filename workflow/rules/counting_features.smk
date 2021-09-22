# count features
rule counting_features:
    input:
        map_file = mapped_reads_dir / "{sample}_mapped_sorted.bam",
        map_index = mapped_reads_dir / "{sample}_mapped_sorted.bam.bai",
        gff = concat_gff_mod_file,
    output:
        feature_count_dir / "{sample}.tsv",
    resources:
        mem_mb=100000,
    conda:
        "../envs/htseq-count.yaml"
    shell:
        "htseq-count --idattr=ID -a 10 -s reverse --nonunique all "
        "{input.map_file} {input.gff} > {output}"
