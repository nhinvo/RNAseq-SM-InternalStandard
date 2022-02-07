# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = output_path_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = output_path_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        ref = output_path_dict["concat_genome"]["concat_genome_file"],
        indexing = output_path_dict["concat_genome"]["concat_genome_done"]
    output:
        sam_out = temp(output_path_dict["mapped_reads"] / "{sample}_mapped.sam"),
    resources: 
        mem_mb=100000,
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output.sam_out}
        """

rule index_genome:
    input:
        ref = output_path_dict["concat_genome"]["concat_genome_file"]
    output:
        touch(output_path_dict["concat_genome"]["concat_genome_done"]),
    resources:
        mem_mb=10000,
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build {input.ref} {input.ref}"

