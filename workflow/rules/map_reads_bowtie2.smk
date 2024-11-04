# combine all reference genome - for read mapping 
rule combine_genome_references:
    input:
        [d['genome'] for d in config["input"]["reference genomes"].values()]
    output:
        scratch_dict["concat_genome"]["concat_genome_file"]
    shell:
        "cat {input} > {output}; sleep 1"

# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = scratch_dict["concat_genome"]["concat_genome_file"],
        indexing = scratch_dict["concat_genome"]["concat_genome_done"]
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output} -p {resources.tasks}"

rule index_genome:
    input:
        ref = scratch_dict["concat_genome"]["concat_genome_file"]
    output:
        touch(scratch_dict["concat_genome"]["concat_genome_done"]),
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --threads {resources.tasks} {input.ref} {input.ref}"

