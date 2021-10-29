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

results_dir = Path(snakemake.input.r_dir)
htseq2dir = Path(snakemake.input.sample_counts[0]).parent
raw_gff_dir = Path(snakemake.input.raw_gff_dir)
condition_table_path = Path(snakemake.input.condition_table_path)
raw_reads_dir = Path(snakemake.input.raw_reads)

print("results_dir: {}".format(results_dir))
print("htseq2dir: {}".format(htseq2dir))
print("raw_gff_dir: {}".format(raw_gff_dir))
print("condition_table_path: {}".format(condition_table_path))

main(results_dir, htseq2dir, raw_gff_dir, condition_table_path, raw_reads_dir, feature_types_to_keep)