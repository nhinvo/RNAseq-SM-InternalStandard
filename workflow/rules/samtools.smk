
rule convert_sam2bam:
    input:
        Path(config["output"]["mapped_reads"]) / "{sample}_mapped.sam",
    output:
        temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        Path(config["output"]["mapped_reads"]) / "{sample}_mapped.bam",
    output:
        Path(config["output"]["mapped_reads"]) / "{sample}_mapped_sorted.bam",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule index_bam:
    input:
        Path(config["output"]["mapped_reads"]) / "{sample}_mapped_sorted.bam",
    output:
        temp(Path(config["output"]["mapped_reads"]) / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input}"