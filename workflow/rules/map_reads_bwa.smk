rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = scratch_dict["concat_genome"]["concat_genome_file"],
        indexing = scratch_dict["concat_genome"]["concat_genome_done"]
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'map_reads_bwa' / 'map_reads_PE', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'map_reads_bwa' / 'map_reads_PE', wildcards.sample),
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.r1} {input.r2} > {output}"


rule index_genome:
    input:
        ref = scratch_dict["concat_genome"]["concat_genome_file"]
    output:
        touch(scratch_dict["concat_genome"]["concat_genome_done"]),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = mk_out(log_dir / 'map_reads_bwa' / 'index_genome.out'),
        error = mk_err(log_dir / 'map_reads_bwa' / 'index_genome.err'),
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.ref}"