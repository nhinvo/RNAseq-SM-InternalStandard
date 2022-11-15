
rule run_trim_PE:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'forward read'],
        r2 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'reverse read'],
        ref = Path(config["input"]["adapter_file"]),
    output:
        o1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        o2 = scratch_dict["trimmed_reads"] / "{sample}_2_trimmed.fastq.gz",
    resources:
        partition = 'sched_mit_chisholm',
        mem = '250G',
        ntasks = 20,
        time = '0-12', 
        output = lambda wildcards: mk_out(log_dir / 'run_trim' / 'run_trim_PE', wildcards.sample),
        error = lambda wildcards: mk_err(log_dir / 'run_trim' / 'run_trim_PE', wildcards.sample),
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh threads={resources.ntasks} "
        "in1={input.r1} in2={input.r2} "
        "out1={output.o1} out2={output.o2} "
        "minlen=25 qtrim=rl trimq=10 "
        "ref={input.ref} ktrim=r k=23 mink=11 hdist=1"

