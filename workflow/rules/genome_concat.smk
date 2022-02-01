# concat all genomes within genome directory (necessary for mapping step)
rule genome_concat:
    input:
        config["input"]["genome_refs"],
    output:
        output_path_dict["concat_genome"]["concat_genome_file"],
    resources:
        mem_mb=100000,
    shell:
        "cat {input}/*.fna > {output}; sleep 1"