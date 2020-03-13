import sys
import os

infile = open(sys.argv[1])
outfile1 = open("newRead_1.fastq","w")
outfile2 = open("newRead_2.fastq","w")

while True:
    header1 = infile.readline().rstrip()
    if not header1:
        break
    seq1 = infile.readline().rstrip()
    infile.readline().rstrip()
    qual1 = infile.readline().rstrip()

    header2 = infile.readline().rstrip()
    seq2 = infile.readline().rstrip()
    infile.readline().rstrip()
    qual2 = infile.readline().rstrip()

    if not (header1.split("/"))[0] == (header2.split("/"))[0]:
        print("The reads are not in the same order. Now exit")
        print(header1)
        print(header2)
        exit()

    outfile1.write(header1+"\n"+seq1+"\n+\n"+qual1+"\n")
    outfile2.write(header2+"\n"+seq2+"\n+\n"+qual2+"\n") 

outfile1.close()
outfile2.close() 