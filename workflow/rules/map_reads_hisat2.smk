# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = Path(config["output"]["trimmed_reads"]) / "{sample}_1_trimmed.fastq",
        r2 = Path(config["output"]["trimmed_reads"]) / "{sample}_2_trimmed.fastq",
        index_dir = Path(config["output"]["genome_index_parent"]) / "{sample}" / "concat_genome", 
        index = Path(config["output"]["done_files"]) / "build_hisat2_index.{sample}.done"
    output:
        temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    shell:
        # HISAT2 Mapping:
        """
        hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output}
        """

rule map_reads_SE:
    input:
        se_reads = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.fastq",
        index_dir = Path(config["output"]["genome_index_parent"]) / "{sample}" / "concat_genome", 
        index = Path(config["output"]["done_files"]) / "build_hisat2_index.{sample}.done"
    output:
        temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    shell:
        # HISAT2 Mapping:
        """
        hisat2 -x {input.index_dir} -U {input.se_reads} > {output}
        """

rule make_sample_index_dir:
    output:
        directory(Path(config["output"]["genome_index_parent"]) / "{sample}" / "concat_genome")
    shell:
        "mkdir -p {output}; sleep 1"

rule build_hisat2_index:
    input:
        out_dir = Path(config["output"]["genome_index_parent"]) / "{sample}" / "concat_genome",
        ref = Path(config["output"]["concat_genome"]["concat_genome_file"]),
    output:
        touch(Path(config["output"]["done_files"]) / "build_hisat2_index.{sample}.done")
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    shell:
        "hisat2-build {input.ref} {input.out_dir}"