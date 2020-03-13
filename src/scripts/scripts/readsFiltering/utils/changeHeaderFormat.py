import sys
import os

read1 = sys.argv[1]
read2 = sys.argv[2]

infile1 = open(read1)
infile2 = open(read2)

outfile1 = open("temp1","w")
outfile2 = open("temp2","w")


while True:
    header1 = infile1.readline().rstrip()
    if not header1:
        break
    fields = header1.split(" ")
    newHeader1 = fields[0]+"/1"
    seq1 = infile1.readline().rstrip()
    infile1.readline().rstrip()
    qual1 = infile1.readline().rstrip()

    header2 = infile2.readline().rstrip()
    if not header2:
        break
    fields = header2.split(" ")
    newHeader2 = fields[0]+"/2"
    seq2 = infile2.readline().rstrip()
    infile2.readline().rstrip()
    qual2 = infile2.readline().rstrip()

    outfile1.write(newHeader1+"\n"+seq1+"\n+\n"+qual1+"\n")
    outfile2.write(newHeader2+"\n"+seq2+"\n+\n"+qual2+"\n")

outfile1.close()
outfile2.close()

os.system("mv temp1 "+read1)
os.system("mv temp2 "+read2)
