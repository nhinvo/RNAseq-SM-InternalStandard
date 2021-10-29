from post_htseq2_parsing import *
from pathlib import Path

feature_types_to_keep=[
'CDS',
'rRNA',
'tRNA',
# 'RNase_P_RNA',
# 'pseudogene',
# 'region',
# 'riboswitch',
# 'tmRNA',
]

global_dir = Path("/nobackup1/kve/2021_Sar11ProProject/data")

input_dir = global_dir / "input_data"
output_dir = global_dir / "output_data"
results_dir = global_dir / "results"

raw_reads_dir = input_dir / "raw_reads"
genome_refs_dir = input_dir / "culture_genome_refs"
gff_refs_dir = input_dir / "culture_genome_annotations"
adapter_file = input_dir / "adapters" / "all_illumina_adapters.fa"
condition_table_file = input_dir / "exp_sample_conditions" / "Xinda_Sample_list_GC.tsv"

done_file_dir = output_dir / "done_files"
trimmed_reads_dir = output_dir / "trimmed_reads"
concat_gff_file = output_dir / "concat_gff" / "concat_gff.gff"
concat_gff_mod_file = output_dir / "concat_gff" / "concat_gff_mod.gff"
concat_genome_file = output_dir / "concat_genome" / "concat_genome.fna"
mapped_reads_dir = output_dir / "mapped_reads"
feature_count_dir = output_dir / "HTseq"
genome_index_parent_dir = output_dir / "genome_index"


main(results_dir, feature_count_dir, gff_refs_dir, condition_table_file, raw_reads_dir, feature_types_to_keep)