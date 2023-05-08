
rule convert_sam2bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {resources.tasks} -S -b {input} > {output}"

rule sort_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {resources.tasks} {input} -o {output}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {resources.tasks} -b {input}"

