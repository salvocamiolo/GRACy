import sys
import os
from Bio import SeqIO


genome2test = sys.argv[1]
refGenome = sys.argv[2]

numSnps = 0

os.system("cat "+genome2test+" "+refGenome+" > mafft_input.fasta")
os.system("/home/gracy/Desktop/GRACy_easyinstall/src/conda/bin/mafft mafft_input.fasta > mafft_output.fasta")

numSeq = 0
for seq_record in SeqIO.parse("mafft_output.fasta","fasta"):
    if numSeq == 0:
        sequence2test = str(seq_record.seq)
    else:
        refSequence = str(seq_record.seq)
    numSeq+=1

for a in range(len(sequence2test)):
    if not sequence2test[a] == 'n' and not sequence2test[a]=='-' and not refSequence[a] == '-':
        if not sequence2test[a] == refSequence[a]:
            print(sequence2test[a],refSequence[a],a)
            numSnps += 1

outfile = open("alignmentAssemblyTest_result.txt","w")
outfile.write("SeqLength\tRefLength\tNs_inSeq\tMissing_inSeq\tMissin_inRef\tNum_SNPs\n")
Ns_inSeq = sequence2test.count('n')+sequence2test.count('N')
Missing_inSeq = sequence2test.count('-')
Missing_inRef = refSequence.count('-') #- sequence2test.count('n') - sequence2test.count('N')
outfile.write(str(len(sequence2test))+"\t"+str(len(refSequence))+"\t"+str(Ns_inSeq)+"\t"+str(Missing_inSeq)+"\t"+str(Missing_inRef)+"\t"+str(numSnps)+"\n")
