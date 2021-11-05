rule unzip_gz:
    input:
        Path(config["output"]["trimmed_reads"]) / "{anything}.fastq.gz",
    output:
        Path(config["output"]["trimmed_reads"]) / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"

rule unzip_bz2:
    input:
        Path(config["output"]["trimmed_reads"]) / "{anything}.fastq.bz2",
    output:
        Path(config["output"]["trimmed_reads"]) / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "bzip2 -d {input}"
        