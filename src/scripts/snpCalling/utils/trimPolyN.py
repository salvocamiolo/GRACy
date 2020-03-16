import os
import sys

read1 = sys.argv[1]
read2 = sys.argv[2]

infile1 = open(read1)
infile2 = open(read2)

outfile1 = open( "temp764289_1.fastq","w")
outfile2 = open( "temp764289_2.fastq","w")




while True:
    header1 = infile1.readline().rstrip()
    if not header1:
        break
    seq1 = infile1.readline().rstrip()
    infile1.readline().rstrip()
    qual1 = infile1.readline().rstrip()

    header2 = infile2.readline().rstrip()
    if not header2:
        break
    seq2 = infile2.readline().rstrip()
    infile2.readline().rstrip()
    qual2 = infile2.readline().rstrip()

    #trim read1
    a = len(seq1)
    while True:
        if seq1[a-40:a].count("G")>30 or seq1[a-40:a].count("C")>30:
            a = a -1
        else:
            newSeq1 = seq1[:a-40]
            newQual1 = qual1[:a-40]
            break

    a = len(seq2)
    while True:
        if seq2[a-40:a].count("G")>30 or seq2[a-40:a].count("C")>30:
            a = a -1
        else:
            newSeq2 = seq2[:a-40]
            newQual2 = qual2[:a-40]
            break

    if len(newSeq1) >= 80 and len(newSeq2)>=80:
        outfile1.write(header1+"\n"+newSeq1+"\n+\n"+newQual1+"\n")
        outfile2.write(header2+"\n"+newSeq2+"\n+\n"+newQual2+"\n")

outfile1.close()
outfile2.close()

os.system("mv temp764289_1.fastq "+read1)
os.system("mv temp764289_2.fastq "+read2)
