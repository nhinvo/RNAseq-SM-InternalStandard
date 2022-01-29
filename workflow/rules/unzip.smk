rule unzip_gz:
    input:
        path_dict["trimmed_reads"] / "{anything}.fastq.gz",
    output:
        path_dict["trimmed_reads"] / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"

# rule unzip_bz2:
#     input:
#         path_dict["trimmed_reads"]) / "{anything}.fastq.bz2",
#     output:
#         path_dict["trimmed_reads"]) / "{anything}.fastq",
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/gzip.yaml"
#     shell:
#         "bzip2 -d {input}"
        