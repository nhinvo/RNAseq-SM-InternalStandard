
import gffpandas.gffpandas as gffpd

print(snakemake.input)
print(snakemake.output)

feature_types_to_keep = snakemake.config["feature_types_to_keep"]

annotation = gffpd.read_gff3(str(snakemake.input))

#remove non-desired feature types
annotation.df = annotation.df[annotation.df.type.isin(feature_types_to_keep)]

#create just one gff with feature replaced as exon purely for htseq counting. Then use the unmodifed gff in the post-htseq-parsing.py
annotation.df.iloc[:,2] = 'exon'

annotation.to_gff3(str(snakemake.output))