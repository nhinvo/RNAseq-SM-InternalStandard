# count features
rule counting_features:
    input:
        map_file = mapped_reads_dir / "{sample}_mapped_sorted.bam",
        map_index = mapped_reads_dir / "{sample}_mapped_sorted.bam.bai",
        gff = concat_gff_mod_file,
    output:
        feature_count_dir / "{sample}.tsv",
    resources:
        mem_mb=100000,
    conda:
        "../envs/htseq-count.yaml"
    shell:
        "htseq-count --idattr=ID -a 10 -s reverse --nonunique all "
        "{input.map_file} {input.gff} > {output}"

rule counting_reads:
    input:
        raw_reads_dir
    output:
        library_count_dir / "library_len.tsv"
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