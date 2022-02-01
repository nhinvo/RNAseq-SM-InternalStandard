rule unzip_gz:
    input:
        output_path_dict["trimmed_reads"] / "{anything}.fastq.gz",
    output:
        output_path_dict["trimmed_reads"] / "{anything}.fastq",
    resources:
        mem_mb=100000,
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -d {input}"

# rule unzip_bz2:
#     input:
#         output_path_dict["trimmed_reads"]) / "{anything}.fastq.bz2",
#     output:
#         output_path_dict["trimmed_reads"]) / "{anything}.fastq",
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/gzip.yaml"
#     shell:
#         "bzip2 -d {input}"
        