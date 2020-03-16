
import sys
from Bio import SeqIO
import os

for seq_record in SeqIO.parse(sys.argv[1],"fasta"):
    reference = str(seq_record.seq)


referenceList = []
for a in range(len(reference)):
    referenceList.append(reference[a])


infile = open(sys.argv[2])
outfile = open("Nranges.txt","w")

while True:
    line=infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")

    if int(fields[3])<3 or fields[2]=='N':
        rangeStart = int(fields[1])
        while int(fields[3]) <3 or fields[2]=='N':
            line=infile.readline().rstrip()
            if not line:
                break
            fields = line.split()
        rangeEnd = int(fields[1])
        if rangeStart > 3000 and rangeEnd<(len(reference)-3000):
            outfile.write(str(rangeStart)+"\t"+str(rangeEnd)+"\n")

            if rangeStart>500 and rangeEnd < (len(reference)-500):
                for a in range(rangeStart-500, rangeEnd+500):
                    referenceList[a] = 'N'
    

newReference = "".join(referenceList)
outfile2 = open(sys.argv[1]+"_N.fasta","w")
outfile2.write(">finalScaffold\n"+newReference+"\n")
outfile2.close()

infile.close()

infile = open(sys.argv[2])
line = infile.readline().rstrip()
fields = line.split("\t")
position = int(fields[1])

while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if int(fields[1]) > position + 1:
        outfile.write(str(position)+"\t"+fields[1]+"\n")
        position = int(fields[1])
    else:
        position += 1
        






outfile.close()
os.system("sort -k1n,1n Nranges.txt >temp;mv temp Nranges.txt")
#smooth low coverge Ranges
infile = open("Nranges.txt")
outfile = open("Nranges_smooth.txt","w")

line = infile.readline().rstrip()
if not line:
    exit()
fields = line.split("\t")

value_1a = int(fields[0])
value_1b = int(fields[1])

line = infile.readline().rstrip()
if not line:
    print("One gap to fill")
    outfile.write(str(value_1a)+"\t"+str(value_1b)+"\n")
    exit()
else:
    infile.seek(0)

while True:
    line = infile.readline().rstrip()
    if not line:
        outfile.write(str(value_2a)+"\t"+str(value_2b)+"\n")
        break
    fields = line.split("\t")
    value_2a = int(fields[0])
    value_2b = int(fields[1])
    if (value_2a - value_1b ) < 1000:
        value_1b = value_2b
    else:
        outfile.write(str(value_1a)+"\t"+str(value_1b)+"\n")
        value_1a = value_2a
        value_1b = value_2b







