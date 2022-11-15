
rule convert_sam2bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped.bam"),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'samtools' / 'convert_sam2bam', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'samtools' / 'convert_sam2bam', wildcards.sample),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {resources.ntasks} -S -b {input} > {output}"

rule sort_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'samtools' / 'sort_bam', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'samtools' / 'sort_bam', wildcards.sample),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {resources.ntasks} {input} -o {output}"

rule index_bam:
    input:
        scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(scratch_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'samtools' / 'index_bam', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'samtools' / 'index_bam', wildcards.sample),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {resources.ntasks} -b {input}"

