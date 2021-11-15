# # map sample reads to concatenated genome
# rule map_reads_PE:
#     input:
#         r1 = Path(config["output"]["trimmed_reads"]) / "{sample}_1_trimmed.fastq",
#         r2 = Path(config["output"]["trimmed_reads"]) / "{sample}_2_trimmed.fastq",
#         index_dir = Path(config["output"]["genome_index_parent"]) / "{sample}" / "concat_genome", 
#         index = Path(config["output"]["done_files"]) / "build_hisat2_index.{sample}.done"
#     output:
#         temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
#     resources: 
#         mem_mb=100000,
#     conda:
#         "../envs/hisat2.yaml"
#     shell:
#         # HISAT2 Mapping:
#         """
#         hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output}
#         """

# rule map_reads_SE_greater_70bp:
#     input:
#         se_reads = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.fastq",
#         ref = Path(config["output"]["concat_genome"]["concat_genome_file"]),
#     output:
#         temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/bwa.yaml"
#     shell:
#         # BWA Mapping:
#         """
#         bwa index {input.ref}
#         bwa mem {input.ref} {input.se_reads} > {output}
#         """

rule index_genome:
    input:
        ref = Path(config["output"]["concat_genome"]["concat_genome_file"])
    output:
        touch(Path(config["output"]["concat_genome"]["concat_genome_done"])),
    resources:
        mem_mb=10000,
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.ref}"


rule map_reads_SE_less_70bp:
    input:
        se_reads = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.fastq",
        ref = Path(config["output"]["concat_genome"]["concat_genome_file"]),
        indexing = Path(config["output"]["concat_genome"]["concat_genome_done"])
    output:
        sam = temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
        index = temp(Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.sai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/bwa.yaml"
    shell:
        # BWA Mapping:
        """
        bwa aln {input.ref} {input.se_reads} > {output.index}
        bwa samse {input.ref} {output.index} {input.se_reads} > {output.sam}
        """
