# map sample reads to concatenated genome
rule map_reads_PE:
    input:
        r1 = output_path_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq",
        r2 = output_path_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq",
        index_dir = output_path_dict["genome_index_parent"] / "{sample}" / "concat_genome", 
        index = output_path_dict["done_files"] / "build_hisat2_index.{sample}.done"
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/hisat2/map_reads_PE/{sample}.log"
    shell:
        # HISAT2 Mapping:
        "hisat2 -x {input.index_dir} -1 {input.r1} -2 {input.r2} > {output} 2> {log}"
        

rule map_reads_SE:
    input:
        se_reads = output_path_dict["trimmed_reads"] / "{sample}_trimmed.fastq",
        index_dir = output_path_dict["genome_index_parent"] / "{sample}" / "concat_genome", 
        index = output_path_dict["done_files"] / "build_hisat2_index.{sample}.done"
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.sam"),  # temp() temporarily keeps SAMs until converted
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
        directory(output_path_dict["genome_index_parent"] / "{sample}" / "concat_genome")
    shell:
        "mkdir -p {output}; sleep 1"

rule build_hisat2_index:
    input:
        out_dir = output_path_dict["genome_index_parent"] / "{sample}" / "concat_genome",
        ref = output_path_dict["concat_genome"]["concat_genome_file"],
    output:
        touch(output_path_dict["done_files"] / "build_hisat2_index.{sample}.done")
    resources:
        mem_mb=100000,
    conda:
        "../envs/hisat2.yaml"
    log: 
        "logs/hisat2/build_hisat2_index.log"
    shell:
        "hisat2-build {input.ref} {input.out_dir} &> {log}"
        