
rule convert_sam2bam:
    input:
        path_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(path_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        path_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule index_bam:
    input:
        path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input}"


rule bam_coverage:
    input:
        path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        path_dict["coverage_positions"] / "{sample}_coverage_depth.tsv",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools depth -a -H {input} > {output}"