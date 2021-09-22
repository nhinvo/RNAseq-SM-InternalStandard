from pathlib import Path
from Bio import SeqIO
import pandas as pd
import gzip

def main(raw_reads_dir):
    
    raw_read_library_counts = {}
    
    for i, path in enumerate(sorted(raw_reads_dir.iterdir())):
        sample = path.stem
        counter = 0
        print(sample)
        # with gzip.open(path, "rt") as handle:
        #     first = next(SeqIO.parse(handle, "fastq"))
        #     print(first)
        #     second = next(SeqIO.parse(handle, "fastq"))
        #     print(second)
        # # 
        with gzip.open(path, "rt") as handle:
            counter = len(list(SeqIO.parse(handle, "fastq")))

        print(counter)
        print()
        raw_read_library_counts[sample] = counter
        
    library = pd.DataFrame.from_dict(raw_read_library_counts, orient='index')
    print(library)

global_dir = Path("/nobackup1/kve/2021_Sar11ProProject/data")
input_dir = global_dir / "input_data"
raw_reads_dir = input_dir / "raw_reads"

main(raw_reads_dir)