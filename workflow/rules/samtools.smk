
rule convert_sam2bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/convert_sam2bam/{sample}.log"
    shell:
        "samtools view -S -b {input} > {output} 2> {log}"

rule sort_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/sort_bam/{sample}.log"
    shell:
        "samtools sort {input} -o {output} &> {log}"

rule index_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index_bam/{sample}.log"
    shell:
        "samtools index -b {input} &> {log}"


rule bam_coverage:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        output_path_dict["coverage_positions"] / "{sample}_coverage_depth.tsv",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/samtools/index_bam/{sample}.log"
    shell:
        "samtools depth -a -H {input} > {output} 2> {log}"