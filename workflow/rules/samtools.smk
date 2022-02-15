
rule convert_sam2bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -S -b {input} > {output}"

rule sort_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -n {input} -o {output}"

rule index_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -b {input}"


rule bam_coverage:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        output_path_dict["coverage_positions"] / "{sample}_coverage_depth.tsv",
    resources:
        mem_mb=100000,
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools depth -a -H {input} > {output}"