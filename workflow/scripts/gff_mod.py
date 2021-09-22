
import gffpandas.gffpandas as gffpd

print(snakemake.input)
print(snakemake.output)

annotation = gffpd.read_gff3(str(snakemake.input))

#create just one gff with feature replaced as exon purely for htseq counting. Then use the unmodifed gff in the post-htseq-parsing.py
annotation.df.iloc[:,2] = 'exon'

annotation.to_gff3(str(snakemake.output))