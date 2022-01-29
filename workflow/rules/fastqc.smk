
rule run_fastqc:
    input:
        path_dict["trimmed_reads"] / "{anything}_trimmed.fastq"
    output:
        html = path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input}"

rule move_fastqc_data:
    input:
        html = path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    output:
        html = path_dict["fastqc"] / "{anything}_trimmed_fastqc.html",
        zipp = path_dict["fastqc"] / "{anything}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
