
rule run_trim_PE:
    input:
        r1 = Path(config["input"]["raw_reads"]) / "2021_Sar11Pro_RNAseq_Project_{sample}_1_sequence.fastq.gz",
        r2 = Path(config["input"]["raw_reads"]) / "2021_Sar11Pro_RNAseq_Project_{sample}_2_sequence.fastq.gz",
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = Path(config["output"]["trimmed_reads"]) / "{sample}_1_trimmed.fastq.gz",
        o2 = Path(config["output"]["trimmed_reads"]) / "{sample}_2_trimmed.fastq.gz",
    resources:
        mem_mb=100000,
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"



# rule run_trim_SE:
#     input:
#         r = Path(config["input"]["raw_reads"]) / "2019_dark_adapted_transcriptome_sequencing_{sample}.fastq.bz2",
#         ref = Path(config["input"]["adapter_file"]),
#     output:
#         o = Path(config["output"]["trimmed_reads"]) / "{sample}_trimmed.fastq.gz",
#     resources:
#         mem_mb=100000,
#     conda:
#         "../envs/bbtools.yaml"
#     shell:
#         "bbduk.sh "
#         "in={input.r} "
#         "out={output.o} "
#         "minlen=25 qtrim=rl trimq=10 "
#         "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"