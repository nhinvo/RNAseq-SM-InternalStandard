RAW_READS_DIR="/nobackup1/kve/2021_Sar11ProProject/data/input_data/raw_reads"
echo -e "read_file\tnum_reads" > "library_len.tsv"
for READ in ${RAW_READS_DIR}/*
do
    FILE_BASENAME=$(basename $READ)
    FILE_COUNT=$(zcat $READ | echo $((`wc -l`/4)))
    echo -e "$FILE_BASENAME\t$FILE_COUNT"  >> "library_len.tsv"
    echo $FILE_BASENAME
done