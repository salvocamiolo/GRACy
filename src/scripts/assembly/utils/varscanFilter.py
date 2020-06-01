import os,sys
from Bio import SeqIO
from Bio import Seq
import argparse
parser = argparse.ArgumentParser(description="Filter varscan output in vcf format by only keeping variants with the highest number of reads")
parser.add_argument("-o","--output",required=True,help="The varscn output file")
parser.add_argument("-i","--input",required=True,help="The varscan input file")
parser.add_argument("-g","--validateIndel",required=False,help="Validate indels by searching the flanking motif in the reads")
parser.add_argument("-1","--read1",required=False,help="First fastq file in reads, needed for indel validation")
parser.add_argument("-2","--read2",required=False,help="Second fastq file in reads, needed for indel validation")
parser.add_argument("-r","--reference",required=False,help="Reference fasta file, needed for indel validation")

args = vars(parser.parse_args())
inputFile = args['input']
outputFile = args['output']
reference = args['reference']
if not args['validateIndel']:
    validate = 0
else:
    validate =1

read1 = args['read1']
read2 = args['read2']

kmersInReads = []
if validate == 1:
    print("Reading read1")
    for seq_record in SeqIO.parse(read1,"fastq"):
        sequence = str(seq_record.seq)
        for a in range(0,len(sequence)-50,+50):
            kmersInReads.append(sequence[a:a+50])
    
    print("Reading read2")
    for seq_record in SeqIO.parse(read2,"fastq"):
        sequence = str(seq_record.seq)
        for a in range(0,len(sequence)-50,+50):
            kmersInReads.append(sequence[a:a+50])

    for seq_record in SeqIO.parse(reference,"fasta"):
        refSeq = str(seq_record.seq)
    




infile = open(inputFile)
outfile = open(outputFile,"w")

line = infile.readline().rstrip()
outfile.write(line+"\n")

while not "#CHROM" in line:
    line = infile.readline().rstrip()
    outfile.write(line+"\n")

while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    info = fields[9].split(":")
    if len(fields[3]) >1 or len(fields[4])>1: #an indel
        if validate == 1:
            print("Validating deletion at position %s" %fields[1])
            refAllele = refSeq[int(fields[1])-25:int(fields[1])-2]+"-"+fields[3]+"-"+refSeq[int(fields[1]):int(fields[1])+25]
            altAllele = refSeq[int(fields[1])-25:int(fields[1])-2]+"-"+fields[4]+"-"+refSeq[int(fields[1]):int(fields[1])+25]
            print(refAllele)
            print(altAllele)
            sys.stdin.read(1)
    
    
    if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
        outfile.write(line+"\n")


infile.close()
outfile.close()


    



