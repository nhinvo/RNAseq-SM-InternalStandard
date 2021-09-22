# concat all GFF within GFF directory (necessary for counting step)
rule annotation_concat:
    input:
        gff_refs_dir,
    output:
        concat_gff_file,
    resources:
        mem_mb=100000,
    shell:
        "cat {input}/*.gff > {output}; sleep 1"

rule add_exon_column_to_gff:
    input:
        concat_gff_file
    output:
        concat_gff_mod_file
    conda:
        "../envs/post_htseq2_parsing.yaml"
    resources:
        mem_mb=1000,
    script:
        "../scripts/gff_mod.py"
