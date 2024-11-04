rule convert_sam2bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {resources.cpus_per_task} -S -b {input} > {output}"

rule sort_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {resources.cpus_per_task} {input} -o {output}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {resources.cpus_per_task} -b {input}"

rule mapping_coverage:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        scratch_dict["mapping_coverage"] / "{sample}.cov", 
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools depth {input} > {output}"