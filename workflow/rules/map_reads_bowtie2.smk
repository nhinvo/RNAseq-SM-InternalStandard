# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = output_path_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = output_path_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = output_path_dict["concat_genome"]["concat_genome_file"],
        indexing = output_path_dict["concat_genome"]["concat_genome_done"]
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'map_reads_bowtie2' / 'map_reads_PE', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'map_reads_bowtie2' / 'map_reads_PE', wildcards.sample),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output} -p {resources.ntasks}"

rule index_genome:
    input:
        ref = output_path_dict["concat_genome"]["concat_genome_file"]
    output:
        touch(output_path_dict["concat_genome"]["concat_genome_done"]),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = mk_out(log_dir / 'map_reads_bowtie2' / 'index_genome.out'),
        error = mk_err(log_dir / 'map_reads_bowtie2' / 'index_genome.err'),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --threads {resources.ntasks} {input.ref} {input.ref}"

