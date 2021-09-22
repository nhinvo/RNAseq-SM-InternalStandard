# concat all genomes within genome directory (necessary for mapping step)
rule genome_concat:
    input:
        genome_refs_dir,
    output:
        concat_genome_file,
    resources:
        mem_mb=100000,
    shell:
        "cat {input}/*.fna > {output}; sleep 1"