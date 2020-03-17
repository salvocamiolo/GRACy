import os
import sys
from Bio import SeqIO


vcfFileName = sys.argv[1]
cdsFile = sys.argv[2]
gffFile = sys.argv[3] #In the GUI version this comes from input



#Create codon Table
bases = "tcag"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_tableL = dict(zip(codons, amino_acids))

bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_tableC = dict(zip(codons, amino_acids))

codon_table = {}
codon_table =codon_tableC.copy()
codon_table.update(codon_tableL)
#print(codon_table)



#Parse gff file
infile = open(gffFile)
cdsCoord = {}
geneCoord = {}
strand = {}
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if len(fields) >3:
        if fields[2]=="gene":
            geneName = ((fields[8].split("Name="))[1].split(";"))[0]

            if not geneName in strand:
                strand[geneName] = fields[6]
            if not geneName in cdsCoord:
                cdsCoord[geneName] = []
                strand[geneName] = fields[6]
            line2 = infile.readline().rstrip()
            if not line2:
                break
            fields2 = line2.split("\t")
            if fields2[2] == "CDS":
                cdsCoord[geneName].append((int(fields2[3]),int(fields2[4])) )
            filePos = infile.tell()
            while not fields2[2] == "gene":
                line2 = infile.readline().rstrip()
                if not line2:
                    break
                fields2 = line2.split("\t")
                if fields2[2] == "CDS":
                    cdsCoord[geneName].append((int(fields2[3]),int(fields2[4])) )
            infile.seek(filePos)
            if len(cdsCoord[geneName])==0:
                del cdsCoord[geneName]
            else:
                coordinatesSet = []
                for item in cdsCoord[geneName]:
                    coordinatesSet.append(item[0])
                    coordinatesSet.append(item[1])
                    fiveprimeCoord = min(coordinatesSet)
                    threeprimeCoord = max(coordinatesSet)
            if not geneName in geneCoord and  geneName in cdsCoord:
                geneCoord[geneName] = (fiveprimeCoord,threeprimeCoord)



#for item in geneCoord:
#    print item,geneCoord[item]
#    if item in cdsCoord:
#        print cdsCoord[item]


#Collect cds sequences
cdsSequences = {}
for seq_record in SeqIO.parse(cdsFile,"fasta"):
    locus = str(seq_record.id)
    if not locus in cdsSequences:
        cdsSequences[locus] = str(seq_record.seq)






#Analyze vcf file for each experiment in the list

infile.close()

vcfFileList = open(vcfFileName) #To change in GUI version as user input




#*****************************************************************
#****************+ Start main algorithm  *************************
#*****************************************************************

while True:
    experiment = vcfFileList.readline().rstrip()
    if not experiment:
        break

    infile = open(experiment)
    fileName = (experiment.split("/"))[-1]
    pathFileNameParts = (experiment.split("/"))[:-1]
    pathFilename = "/".join(pathFileNameParts)

    outFileName = experiment+"_snpFreq.txt"
    outfile = open(outFileName,"w")
    outfile.write("Position\tRefBase\tTargetBase\tCoverage\tFrequencyRef\tFrequencyTarget\n")

    outSEFileName = experiment+"_snpEffect.txt"
    outfileSE = open(outSEFileName,"w")
    outfileSE.write("Position\tStrand\tRelativePosition\tFrequency\tCoverage\tRefBase\tTargetBase\tRefCodon\tTargetCodon\tRefAA\tTargetAA\n")



    snpReads = {}

    #Read the header
    while not "#CHROM" in line:
        line = infile.readline().rstrip()
    #Scann the file
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        pos = int(fields[1])

        if not pos in snpReads:
            snpReads[pos] = []

        #Calculate frequency
        freqUnits = ((fields[-1].split("="))[-1]).split(",")
        frequency = ((float(freqUnits[2]) + float(freqUnits[3]))/(float(freqUnits[0]) + float(freqUnits[1])+float(freqUnits[2]) + float(freqUnits[3])))*100
        coverage = float(freqUnits[0]) + float(freqUnits[1])+float(freqUnits[2]) + float(freqUnits[3])
        snpReads[pos].append((fields[3],fields[4],float(freqUnits[0]) + float(freqUnits[1]),float(freqUnits[2]) + float(freqUnits[3]) ))



        for item in geneCoord:
            if pos > geneCoord[item][0] and pos < geneCoord[item][1]:
                exonSNP = 0
                if len(item)>0:
                    print("position",pos,"is in gene",item)
                    if strand[item]=="+":
                        relativePosition = 0
                        for ranges in cdsCoord[item]:
                            if ranges[1] >= pos and ranges[0]<=pos:
                                relativePosition +=  pos -ranges[0] + 1
                                exonSNP = 1
                                break
                            else:
                                relativePosition += ranges[1] -ranges[0] +1
                                print("First exonLenfth",relativePosition)
                        print(relativePosition)

                    else:
                        relativePosition = 0
                        introns = 0
                        for ranges in cdsCoord[item]:
                            if ranges[0] <= pos and ranges[1]>=pos:
                                relativePosition +=  -pos +ranges[1] + 1
                                exonSNP = 1
                                break
                            else:
                                relativePosition += ranges[1] -ranges[0] +1
                                print("First exonLenfth",relativePosition)
                                introns =1
                        print(relativePosition)



                    print("The mutated base is on gene",item)
                    print("The position mutated is at",relativePosition)
                    if exonSNP==1 and len(fields[3])==1 and len(fields[4])==1: #Only SNPs are considered at this stage
                        print("Corresponding to the base",cdsSequences[item][relativePosition-1],"in the region",cdsSequences[item][relativePosition-6:relativePosition+4])
                        if strand[item]=="+":
                            newBase = fields[4]
                            refBase = fields[3]
                            if relativePosition%3 == 0:
                                referenceCodon = cdsSequences[item][relativePosition-3:relativePosition]
                                targetCodon = referenceCodon[:2]+fields[4]
                            if relativePosition%3 == 1:
                                referenceCodon = cdsSequences[item][relativePosition-1:relativePosition+2]
                                targetCodon = fields[4]+referenceCodon[1:]
                            if relativePosition%3 == 2:
                                referenceCodon = cdsSequences[item][relativePosition-2:relativePosition+1]
                                targetCodon = referenceCodon[0]+  fields[4] + referenceCodon[2]
                            print("The snp resides in the codon",referenceCodon)

                            print("Target codon",targetCodon)


                        else:
                            if fields[4] == "A" or fields[4] == "a":
                                newBase = "T"
                            if fields[4] == "T" or fields[4] == "t":
                                newBase = "A"
                            if fields[4] == "G" or fields[4] == "g":
                                newBase = "C"
                            if fields[4] == "C" or fields[4] == "c":
                                newBase = "G"


                            if fields[3] == "A" or fields[3] == "a":
                                refBase = "T"
                            if fields[3] == "T" or fields[3] == "t":
                                refBase = "A"
                            if fields[3] == "G" or fields[3] == "g":
                                refBase = "C"
                            if fields[3] == "C" or fields[3] == "c":
                                refBase = "G"



                            if relativePosition%3 == 0:
                                referenceCodon = cdsSequences[item][relativePosition-3:relativePosition]
                                targetCodon = referenceCodon[:2]+newBase
                            if relativePosition%3 == 1:
                                referenceCodon = cdsSequences[item][relativePosition-1:relativePosition+2]
                                targetCodon = newBase+referenceCodon[1:]
                            if relativePosition%3 == 2:
                                referenceCodon = cdsSequences[item][relativePosition-2:relativePosition+1]
                                targetCodon = referenceCodon[0]+  newBase + referenceCodon[2]
                            print("The snp resides in the codon",referenceCodon)

                            print("Target codon",targetCodon)

                        if targetCodon in codon_table and referenceCodon in codon_table:
                            outfileSE.write(item+"\t"+str(pos)+"\t"+strand[item]+"\t"+str(relativePosition)+"\t"+str(frequency)+"\t"+str(coverage)+"\t"+refBase+"\t"+ newBase +"\t"+referenceCodon+"\t"+targetCodon+"\t"+codon_table[referenceCodon]+"\t"+codon_table[targetCodon]+"\n")
                        else:
                            outfileSE.write(item+"\t"+str(pos)+"\t"+strand[item]+"\t"+str(relativePosition)+"\t"+str(frequency)+"\t"+str(coverage)+"\t"+refBase+"\t"+ newBase +"\t"+referenceCodon+"\t"+targetCodon+"\t"+"N_containing_codon"+"\t"+"N_containing_codon"+"\n")






                    else:
                        print("This SNPs is probbaly into an intron")
                    print("The record in the vcf file is",fields[3],fields[4])


    #Calculate frequency of each snp
    snpFreq = {}
    for position in snpReads:
        if not position in snpFreq:
            snpFreq[position] = []
        refReads = snpReads[position][0][2]
        alternativeReads = 0
        for snp in snpReads[position]:
            alternativeReads += snp[3]
        allReads = alternativeReads + refReads
        for snp in snpReads[position]:
            snpFreq[position].append((snp[0],snp[1],snp[2]/allReads,  snp[3]/allReads,snp[2]+snp[3]))


    for position in snpFreq:
        for snp in snpFreq[position]:
            outfile.write(str(position)+"\t"+str(snp[0])+"\t"+str(snp[1])+"\t"+str(snp[4])+"\t"+str(float(snp[2])*100)[:5]+"\t"+str(float(snp[3])*100)[:5]+"\n")


    infile.close()
    outfile.close()
