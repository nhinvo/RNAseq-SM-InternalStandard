rule annotation_concat:
    input:
        [d['annotation'] for d in config["input"]["reference genomes"].values()]
    output:
        scratch_dict["concat_gff"]["concat_gff_file"],
    shell:
        "cat {input} > {output}; sleep 1"

rule add_exon_column_to_gff:
    input:
        scratch_dict["concat_gff"]["concat_gff_file"]
    output:
        scratch_dict["concat_gff"]["concat_gff_mod_file"]
    conda:
        "../envs/uncorrected_count_summary.yaml"
    script:
        "../scripts/gff_mod.py"
