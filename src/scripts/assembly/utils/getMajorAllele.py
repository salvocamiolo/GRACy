import sys
import os


inputFile = sys.argv[1]
outputFile = sys.argv[2]

infile = open(inputFile)
outfile = open(outputFile,"w")

vcfEntry = {}
entryFreq = {}


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
    freq = float((((fields[7].split("AF="))[-1]).split(";"))[0])
    dp4 = ((fields[7].split("DP4="))[-1]).split(",")
    freq = (float(dp4[2])+float(dp[3]))/(float(dp4[2])+float(dp[3])+float(dp4[0])+float(dp[1]))


    if freq > 0.5:
        if not fields[1] in entryFreq:
            entryFreq[fields[1]] = freq
            vcfEntry[fields[1]] = line
        if freq > entryFreq[fields[1]]:
            entryFreq[fields[1]] = freq
            vcfEntry[fields[1]] = line

for item in vcfEntry:
    outfile.write(vcfEntry[item]+"\n")


    
