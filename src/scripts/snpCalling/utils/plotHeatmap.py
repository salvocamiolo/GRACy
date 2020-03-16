
import numpy as np
import sys

position2PlotFile = sys.argv[1]
file2AnalyzeFile = sys.argv[2]

pos2plot = open(position2PlotFile)

outfile = open("snpsTable.txt","w")
outfile.write("Snp\t")

#Read the positions to plot
p2p = []
while True:
    line = pos2plot.readline().rstrip()
    if not line:
        break
    p2p.append(int(line))

p2p_sort = sorted(p2p)

heatmapFiles = open(file2AnalyzeFile)
snpPositions = {}
xlab = []
while True:
    line = heatmapFiles.readline().rstrip()
    if not line:
        break
    
    
    if not line in snpPositions:
        snpPositions[line] = {}
    infile = open(line)
    line2 = infile.readline().rstrip()
    while True:
        line2 = infile.readline().rstrip()
        if not line2:
            break
        fields = line2.split("\t")
        if not (fields[0],fields[1],fields[2]) in snpPositions[line]:
            snpPositions[line][(fields[0],fields[1],fields[2])] = float(fields[5])

    infile.close()
for exp in snpPositions:
    outfile.write(exp+"\t")
    xlab.append(exp)

outfile.write("\n")

#Calculate the total number of snps present in the experiments
allSNPs = set()

for experiments in snpPositions:
    print "In the experiment",experiments,"There are",len(snpPositions[experiments]),"snps"
    for snp in snpPositions[experiments]:
        allSNPs.add(snp)

snp2Plot = []
for position in p2p_sort:
    for item in allSNPs:
        if int(item[0]) == position:
            snp2Plot.append(item)


#Create the table to plot
table2Plot = []
ylab=[]
for snp in snp2Plot:
    ylab.append(snp)
    outfile.write(str(snp)+"\t")
    line = []
    for experiments in snpPositions:
        if snp in snpPositions[experiments]:
            line.append(snpPositions[experiments][snp])
            outfile.write(str(snpPositions[experiments][snp])+"\t")
        else:
            line.append(np.nan)
            outfile.write("\t")
    table2Plot.append(line)
    outfile.write("\n")


    






