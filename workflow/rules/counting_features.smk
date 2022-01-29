# count features
rule counting_features:
    input:
        map_file = Path(config["output"]["mapped_reads"]) / "{sample}_mapped_sorted.bam",
        map_index = Path(config["output"]["mapped_reads"]) / "{sample}_mapped_sorted.bam.bai",
        gff = Path(config["output"]["concat_gff"]["concat_gff_mod_file"]),
    output:
        Path(config["output"]["feature_count"]) / "{sample}.tsv",
    resources:
        mem_mb=10000,
    conda:
        "../envs/htseq-count.yaml"
    params:
        mapping_quality_thresh = config["htseq"]["mapping_qual_threshold"]
    shell:
        "htseq-count --idattr=ID -a {params.mapping_quality_thresh} -s reverse --nonunique all "
        "{input.map_file} {input.gff} > {output}"

rule counting_reads:
    input:
        config["input"]["raw_reads"]
    output:
        config["output"]["library_count"]
    resources:
        mem_mb=100000,
    shell:
        """
        echo -e "read_file\tnum_reads" > {output}
        for READ in {input}/*
        do
            FILE_BASENAME=$(basename $READ)
            FILE_COUNT=$(zcat $READ | echo $((`wc -l`/4)))
            echo -e "$FILE_BASENAME\t$FILE_COUNT"  >> {output}
            echo $FILE_BASENAME
        done
        """