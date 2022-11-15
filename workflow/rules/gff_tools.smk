rule annotation_concat:
    input:
        [d['annotation'] for d in config["input"]["reference genomes"].values()]
    output:
        scratch_dict["concat_gff"]["concat_gff_file"],
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = mk_out(log_dir / 'gff_tools' / 'annotation_concat.out'),
        error = mk_err(log_dir / 'gff_tools' / 'annotation_concat.err'),
    shell:
        "cat {input} > {output}; sleep 1"

rule add_exon_column_to_gff:
    input:
        scratch_dict["concat_gff"]["concat_gff_file"]
    output:
        scratch_dict["concat_gff"]["concat_gff_mod_file"]
    conda:
        "../envs/post_count_analysis.yaml"
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = mk_out(log_dir / 'gff_tools' / 'add_exon_column_to_gff.out'),
        error = mk_err(log_dir / 'gff_tools' / 'add_exon_column_to_gff.err'),
    script:
        "../scripts/gff_mod.py"
