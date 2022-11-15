import pandas as pd

feature_types_to_keep = snakemake.config["feature_types_to_count"]

annotation = pd.read_table(snakemake.input[0], header=None, comment='#')

annotation = annotation[annotation[2].isin(feature_types_to_keep)]

annotation[2] = 'exon'

annotation.to_csv(snakemake.output[0], sep='\t', index=False, header=False)