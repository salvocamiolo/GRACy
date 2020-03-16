
import numpy as np

import sys


file2AnalyzeFile = sys.argv[1]
position2Plot = set()
snpsPerSample = {}

filesToOpen = open(file2AnalyzeFile)

while True:
    sample = filesToOpen.readline().rstrip()
    if not sample:
        break

    infile = open(sample)
    infile.readline()

    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        position2Plot.add((fields[0],fields[1],fields[3],fields[8],fields[9],fields[10],fields[11]))
        if not sample in snpsPerSample:
            snpsPerSample[sample] = {}
        if not (fields[0],fields[1],fields[3],fields[8],fields[9],fields[10],fields[11]) in snpsPerSample[sample]:
            snpsPerSample[sample][(fields[0],fields[1],fields[3],fields[8],fields[9],fields[10],fields[11])] = fields[4]
        

outfile = open(file2AnalyzeFile+"_snpTable.txt","w")
outfile.write("Sample\tGene\tPosGenome\tPosGene\tRefCodon\tAltCodon\t")
for a in snpsPerSample:
    outfile.write(str(a)+"\t")
outfile.write("\n")

for item in position2Plot:
    for item2 in item:
        outfile.write(item2+"\t")
    for item3 in snpsPerSample:
        if item in snpsPerSample[item3]:
            outfile.write(snpsPerSample[item3][item]+"\t")
        else:
            outfile.write("-\t")
    outfile.write("\n")
