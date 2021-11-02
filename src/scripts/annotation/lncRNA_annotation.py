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

for seq_record in SeqIO.parse(inputGenome,"fasta"):
    genomeSeq = str(seq_record.seq)
    genomeID = str(seq_record.id)


outfile1 = open(outputFile+"_lncrna_seqs.fasta","w")
outfile2 = open(outputFile+"_lncrna_gff3.gff","w")

for seq_record in SeqIO.parse(lncRNAfile,"fasta"):
    sequence = str(seq_record.seq)
    seqID = str(seq_record.id)
    seqDescription = str(seq_record.description)
    tempFile = open("tempFasta.fasta","w")
    tempFile.write(">"+seqID+"\n"+sequence)
    tempFile.close()


    os.system(condaFolder+"/bin/makeblastdb -dbtype nucl -in "+inputGenome)
    os.system(condaFolder+"/bin/blastn -query tempFasta.fasta -db "+inputGenome+" -outfmt 6 -task blastn -dust no -soft_masking false > preOutputBlast.txt")
    os.system("head -20 preOutputBlast.txt > preOutputBlast2.txt")
    os.system("sort -k 4rn,4rn -k12rn,12rn preOutputBlast2.txt > outputBlast.txt")
    os.system("rm preOutputBlast.txt preOutputBlast2.txt -f ")
    blastFile = open("outputBlast.txt")
    line = blastFile.readline().rstrip()
    if not line:
        pass
    else:
        fields = line.split()
        coverage = (float(abs(int(fields[9])-int(fields[8])))/float(len(sequence)))*100
        if int(fields[8]) > int(fields[9]):
            strand = "-"
        else:
            strand = "+"
        print(coverage,str(int(fields[9])-int(fields[8])))
        if coverage > 80:
            if strand == "+":
                outfile1.write(">"+seqID+" "+strand+"\n"+genomeSeq[int(fields[8])-1:int(fields[9])]+"\n")
                outfile2.write(genomeID+"\tGRACy\tgene\t"+fields[8]+"\t"+fields[9]+"\t.\t+\t.\tID="+seqID+"_gene;Name="+seqID+";Product="+seqID+"\n")
                outfile2.write(genomeID+"\tGRACy\tlncRNA\t"+fields[8]+"\t"+fields[9]+"\t.\t+\t.\tID="+seqID+"_mRNA;Parent="+seqID+"_gene;Name="+seqID+".1;Product="+seqID+"\n")
            else:
                print(fields[9],fields[8])
                outfile1.write(">"+seqID+" "+strand+"\n"+genomeSeq[int(fields[9])-1:int(fields[8])]+"\n")
                outfile2.write(genomeID+"\tGRACy\tgene\t"+fields[9]+"\t"+fields[8]+"\t.\t-\t.\tID="+seqID+"_gene;Name="+seqID+";Product="+seqID+"\n")
                outfile2.write(genomeID+"\tGRACy\tlncRNA\t"+fields[9]+"\t"+fields[8]+"\t.\t-\t.\tID="+seqID+"_mRNA;Parent="+seqID+"_gene;Name="+seqID+".1;Product="+seqID+"\n")
        else:
            pass

       

     

    blastFile.close()


outfile1.close()


