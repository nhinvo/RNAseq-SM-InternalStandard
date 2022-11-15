
rule run_fastqc:
    input:
        scratch_dict["trimmed_reads"] / "{anything}_trimmed.fastq"
    output:
        html = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'fastqc' / 'run_fastqc', wildcards.anything),
        error = lambda wildcards: mk_err(log_dir / 'fastqc' / 'run_fastqc', wildcards.anything),
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input}"

rule move_fastqc_data:
    input:
        html = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    output:
        html = scratch_dict["fastqc"] / "{anything}_trimmed_fastqc.html",
        zipp = scratch_dict["fastqc"] / "{anything}_trimmed_fastqc.zip",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'fastqc' / 'move_fastqc_data', wildcards.anything),
        error = lambda wildcards: mk_err(log_dir / 'fastqc' / 'move_fastqc_data', wildcards.anything),
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
