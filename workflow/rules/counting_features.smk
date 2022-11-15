# count features
rule counting_features:
    input:
        map_file = scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
        map_index = scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai",
        gff = scratch_dict["concat_gff"]["concat_gff_mod_file"],
    output:
        scratch_dict["feature_count"] / "{sample}.tsv",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '12G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'counting_features' / 'counting_features', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'counting_features' / 'counting_features', wildcards.sample),
    conda:
        "../envs/htseq-count.yaml"
    params:
        mapping_quality_thresh = config["htseq"]["mapping_qual_threshold"]
    shell:
        "htseq-count --idattr=ID -a {params.mapping_quality_thresh} -s reverse -t exon -r pos --nonunique all "
        "{input.map_file} {input.gff} > {output}"
        