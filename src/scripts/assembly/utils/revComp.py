
from Bio import SeqIO
from Bio import Seq

import sys 

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    print(">"+str(seq_record.id)+"\n"+Seq.reverse_complement(seq_record.seq)+"\n")

