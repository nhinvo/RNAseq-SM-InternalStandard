
rule convert_sam2bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.sam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped.bam"),
    benchmark:
        benchmark_dir / 'samtools' / 'convert_sam2bam' / '{sample}.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = str(log_dir / 'samtools' / 'convert_sam2bam' / '{sample}.out'),
        error = str(log_dir / 'samtools' / 'convert_sam2bam' / '{sample}.err'),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -@ {resources.ntasks} -S -b {input} > {output}"

rule sort_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped.bam",
    output:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    benchmark:
        benchmark_dir / 'samtools' / 'sort_bam' / '{sample}.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = str(log_dir / 'samtools' / 'sort_bam' / '{sample}.out'),
        error = str(log_dir / 'samtools' / 'sort_bam' / '{sample}.err'),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {resources.ntasks} {input} -o {output}"

rule index_bam:
    input:
        output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
    output:
        temp(output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam.bai"),
    benchmark:
        benchmark_dir / 'samtools' / 'index_bam' / '{sample}.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = str(log_dir / 'samtools' / 'index_bam' / '{sample}.out'),
        error = str(log_dir / 'samtools' / 'index_bam' / '{sample}.err'),
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {resources.ntasks} -b {input}"


# rule bam_coverage:
#     input:
#         output_path_dict["mapped_reads"] / "{sample}_mapped_sorted.bam",
#     output:
#         output_path_dict["coverage_positions"] / "{sample}_coverage_depth.tsv",
#     benchmark:
#         "benchmark/samtools/bam_coverage/{sample}.benchmark"
#     resources:
#         partition = 'sched_mit_chisholm',
#         mem = '250G',
#         ntasks = 20,
#         time = '0-12', 
#         output = str(log_dir / 'samtools/bam_coverage/{sample}.out'),
#         error = str(log_dir / 'samtools/bam_coverage/{sample}.err'),
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools depth -@ {resources.ntasks} -a -H {input} > {output}"