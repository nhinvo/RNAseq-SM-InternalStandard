
rule run_fastqc:
    input:
        scratch_dict["trimmed_reads"] / "{anything}_trimmed.fastq"
    output:
        html = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.html",
        zipp = scratch_dict["trimmed_reads"] / "{anything}_trimmed_fastqc.zip",
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
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mv {input.html} {output.html}
        mv {input.zipp} {output.zipp}
        """
