import os
import sys
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Annotated the lncRNA in gff files")
parser.add_argument("-f","--fasta",required=True,help="The fasta file with the lncRNA sequences")
parser.add_argument("-i","--inputFile",required=True,help="The fasta file of the genome sequence to annotated")
parser.add_argument("-o","--outputFile",required=True,help="The output file name")
parser.add_argument("-c","--condaFolder",required=True,help="The path to the anaconda/miniconda folder")

args = vars(parser.parse_args())

lncRNAfile = args['fasta']
inputGenome = args['inputFile']
outputFile = args['outputFile']
condaFolder = args['condaFolder']


outfile = open(outputFile,"w")

for seq_record in SeqIO.parse(lncRNAfile,"fasta"):
    sequence = str(seq_record.seq)
    seqID = str(seq_record.id)
    tempFile = open("tempFasta.fasta","w")
    tempFile.write(">"+seqID+"\n"+sequence)
    tempFile.close()


    os.system(condaFolder+"/bin/makeblastdb -dbtype nucl -in "+inputGenome)
    os.system(condaFolder+"/bin/blastn -query tempFasta.fasta -db "+inputGenome+" -outfmt 6 > preOutputBlast.txt")
    os.system("head -20 preOutputBlast.txt > preOutputBlast2.txt")
    os.system("sort -k 3rn,3rn -k12rn,12rn preOutputBlast2.txt > outputBlast.txt")
    os.system("rm preOutputBlast.txt preOutputBlast2.txt -f ")
    blastFile = open("outputBlast.txt")
    print(blastFile.readline())
    blastFile.close()





