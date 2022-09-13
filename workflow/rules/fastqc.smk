
rule run_fastqc:
    input:
        output_path_dict["trimmed_reads"] / "{anything}_trimmed.fastq"
    output:
        html = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc/run_fastqc/{anything}.log"
    shell:
        "fastqc {input}"

rule move_fastqc_data:
    input:
        html = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
    output:
        html = output_path_dict["fastqc"] / "{anything}_trimmed_fastqc.html",
        zipp = output_path_dict["fastqc"] / "{anything}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
