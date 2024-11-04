import pandas as pd
from Bio import SeqIO

internal_seqs_path="/nfs/chisholmlab001/nvo/projects/23_rogier_MIT9301_diel_cycle/DGE-snakemake.01/inputs/ref_genomes/internal_standard_seqs/internal_standard_sequences.txt"
gff_out="/nfs/chisholmlab001/nvo/projects/23_rogier_MIT9301_diel_cycle/DGE-snakemake.01/inputs/ref_genomes/internal_standard_seqs/internal_standard_sequences.gff"

data={'seqname':[],  # name of chromosome
      'source':[], 
      'feature':[],
      'start':[],
      'end':[],
      'score':[],
      'strand':[],
      'frame':[],
      'attribute':[],}
for record in SeqIO.parse(internal_seqs_path, "fasta"):
    data['seqname'].append(record.id) 
    data['source'].append('Gifford 2016')
    data['feature'].append('CDS')
    data['start'].append(1)
    data['end'].append(len(record.seq))
    data['score'].append('.')
    data['strand'].append('+')
    data['frame'].append(0)
    data['attribute'].append('unknown')

df=pd.DataFrame(data)
df.to_csv(gff_out, sep='\t', index=False, header=False)
