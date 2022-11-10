
rule run_fastqc:
    input:
        output_path_dict["trimmed_reads"] / "{anything}_trimmed.fastq"
    output:
        html = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    benchmark:
        benchmark_dir / 'fastqc' / 'run_fastqc' / '{anything}.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = str(log_dir / 'fastqc' / 'run_fastqc' / '{anything}.out'),
        error = str(log_dir / 'fastqc' / 'run_fastqc' / '{anything}.err'),
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input}"

rule move_fastqc_data:
    input:
        html = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    output:
        html = output_path_dict["fastqc"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["fastqc"] / "{anything}_trimmed_fastqc.zip",
    benchmark:
        benchmark_dir / 'fastqc' / 'move_fastqc_data' / '{anything}.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = str(log_dir / 'fastqc' / 'move_fastqc_data' / '{anything}.out'),
        error = str(log_dir / 'fastqc' / 'move_fastqc_data' / '{anything}.err'),
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
