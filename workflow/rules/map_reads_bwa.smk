rule map_reads_PE:
    input:
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = scratch_dict["concat_genome"]["concat_genome_file"],
        indexing = scratch_dict["concat_genome"]["concat_genome_done"]
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.sam"),
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.r1} {input.r2} > {output}"


rule index_genome:
    input:
        ref = scratch_dict["concat_genome"]["concat_genome_file"]
    output:
        touch(scratch_dict["concat_genome"]["concat_genome_done"]),
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.ref}"