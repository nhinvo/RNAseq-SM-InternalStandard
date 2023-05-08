rule unzip_gz:
    input:
        scratch_dict["trimmed_reads"] / "{anything}.fastq.gz",
    output:
        scratch_dict["trimmed_reads"] / "{anything}.fastq",
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"

        