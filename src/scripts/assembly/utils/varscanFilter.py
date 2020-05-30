import os,sys
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description="Filter varscan output in vcf format by only keeping variants with the highest number of reads")
parser.add_argument("-o","--output",required=True,help="The varscn output file")
parser.add_argument("-i","--input",required=True,help="The varscan input file")

args = vars(parser.parse_args())
inputFile = args['input']
outputFile = args['output']



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
    if int(info[5]) > int(info[4]): #the alternate allele has a highe number of reads than the reference allele
        outfile.write(line+"\n")


infile.close()
outfile.close()


    



