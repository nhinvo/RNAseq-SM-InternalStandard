rule unzip_gz:
    input:
        output_path_dict["trimmed_reads"] / "{anything}.fastq.gz",
    output:
        output_path_dict["trimmed_reads"] / "{anything}.fastq",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'unzip' / 'unzip_gz', wildcards.anything),
        error = lambda wildcards: mk_err(log_dir / 'unzip' / 'unzip_gz', wildcards.anything),
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"

        