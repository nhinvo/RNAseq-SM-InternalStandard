rule combine_genome_references:
    input:
        [d['genome'] for d in config["input"]["reference genomes"].values()]
    output:
        output_path_dict["concat_genome"]["concat_genome_file"]
    benchmark:
        benchmark_dir / 'genome_concat' / 'combine_genome_references.benchmark'
    resources:
        partition = 'sched_mit_chisholm',
        mem = '10G',
        ntasks = 1,
        time = '0-12', 
        output = str(log_dir / 'genome_concat' / 'combine_genome_references.out'),
        error = str(log_dir / 'genome_concat' / 'combine_genome_references.err'),
    shell:
        "cat {input} > {output}; sleep 1"