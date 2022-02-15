# count features
rule counting_features:
    input:
        map_file = output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
        map_index = output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai",
        gff = output_path_dict["concat_gff"]["concat_gff_mod_file"],
    output:
        output_path_dict["feature_count"] / "{sample}.tsv",
    resources:
        mem_mb=10000,
    conda:
        "../envs/htseq-count.yaml"
    params:
        mapping_quality_thresh = config["htseq"]["mapping_qual_threshold"]
    shell:
        "htseq-count --idattr=ID -a 10 -s reverse -t exon -r pos --nonunique all "
        "{input.map_file} {input.gff} > {output}"
        