
rule run_fastqc:
    input:
        Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.fastq"
    output:
        html = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed_fastqc.html",
        zipp = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc {input}"

rule move_fastqc_data:
    input:
        html = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed_fastqc.html",
        zipp = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed_fastqc.zip",
    output:
        html = Path(config["output"]["fastqc"]) / "{sample}_trimmed_fastqc.html",
        zipp = Path(config["output"]["fastqc"]) / "{sample}_trimmed_fastqc.zip",
    resources:
        mem_mb=10000,
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
