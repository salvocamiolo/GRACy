from Bio import SeqIO
from Bio.Seq import Seq
import sys

filename = sys.argv[2]+"_"+sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+".txt"
outfile = open(filename,"w")

sequences = {}

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    if not str(seq_record.id) in sequences:
        sequences[str(seq_record.id)] = str(seq_record.seq)
    
sequence = Seq(sequences[sys.argv[2]][int(sys.argv[3])-1:int(sys.argv[4])-1])
if sys.argv[5] == 'r':
    sequence = sequence.reverse_complement()
outfile.write(">"+sys.argv[2]+"_"+sys.argv[3]+"-"+sys.argv[4]+"\n"+str(sequence))
