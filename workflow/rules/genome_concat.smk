rule combine_genome_references:
    input:
        [d['genome'] for d in config["input"]["reference genomes"].values()]
    output:
        scratch_dict["concat_genome"]["concat_genome_file"]
    shell:
        "cat {input} > {output}; sleep 1"