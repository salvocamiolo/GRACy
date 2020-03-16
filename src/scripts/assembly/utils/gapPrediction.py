import os
import sys
import numpy as nm
from Bio import SeqIO


inputFile = sys.argv[1]
read1 = sys.argv[2]
read2 = sys.argv[3]
installationDirectory = sys.argv[4]


for seq_record in SeqIO.parse(inputFile,"fasta"):
    seq2fill = str(seq_record.seq)

os.system(installationDirectory+"src/conda/bin/makeblastdb -dbtype nucl -in hcmv_genomes.fasta")
sequences = {}
for seq_record in SeqIO.parse("hcmv_genomes.fasta","fasta"):
    if str(seq_record.id) not in sequences:
        sequences[str(seq_record.id)] = str(seq_record.seq)

a=-1
while a < len(seq2fill)-1:
    a+=1
    if seq2fill[a]=="N":
        gapStart = a
        fiveprime = seq2fill[a-200:a]
        threeprime = "N"
        while "N" in threeprime:
            a+=1
            threeprime = seq2fill[a:a+200]
        gapEnd = a
        print("\n\nFilling gap",gapStart,gapEnd)
  
        t5 = open("temp5.fasta","w")
        t3 = open("temp3.fasta","w")
        t5.write(">fiveprime\n"+fiveprime+"\n")
        t3.write(">threeprime\n"+threeprime+"\n")
        t5.close()
        t3.close()
        os.system(installationDirectory+"src/conda/bin/blastn -query temp5.fasta -db hcmv_genomes.fasta -outfmt 6 -out temp5_outputBlast.txt >null 2>&1")
        os.system(installationDirectory+"src/conda/bin/blastn -query temp3.fasta -db hcmv_genomes.fasta -outfmt 6 -out temp3_outputBlast.txt >null 2>&1")

        best5hits = {}
        best3hits = {}

        t5 = open("temp5_outputBlast.txt")
        while True:
            line = t5.readline().rstrip()
            if not line:
                break
            fields = line.split("\t")
            if int(fields[7])-int(fields[6]) > 100 and not fields[1] in best5hits and int(fields[8])<int(fields[9]):
                best5hits[fields[1]] = (int(fields[8]),int(fields[9]),int(fields[6]),int(fields[7]))
        
        t3 = open("temp3_outputBlast.txt")
        while True:
            line = t3.readline().rstrip()
            if not line:
                break
            fields = line.split("\t")
            if int(fields[7])-int(fields[6]) > 100 and not fields[1] in best3hits and int(fields[8])<int(fields[9]):
                best3hits[fields[1]] = (int(fields[8]),int(fields[9]),int(fields[6]),int(fields[7]))

        numMatches = 0
        foundSequences = {}
        for item in best5hits:
            if item in best3hits:
                numMatches+=1
                foundSequences[item] = sequences[item][best5hits[item][0]:best3hits[item][1]]
        t5.close()
        t3.close()
        print("Found",numMatches,"sequences featuring similar flanking regions")
        if not numMatches == 0:
            t = open("foundSequences.fasta","w")
            for item in foundSequences:
                t.write(">"+item+"\n"+foundSequences[item]+"\n")
            t.close()
            print("Aligning the reads to the found sequences....")
            os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build foundSequences.fasta found >null 2>&1")
            os.system(installationDirectory+"src/conda/bin/bowtie2 --local -x found -1 "+read1+" -2 "+read2+" -S alignment.sam >null 2>&1")
            print("Converting sam to bam....")
            os.system(installationDirectory+"src/conda/bin/samtools view -bS -h -F 4 alignment.sam >alignment.bam  2>&1")
            print("sorting bam....")
            os.system(installationDirectory+"src/conda/bin/samtools sort -o alignment_sorted.bam alignment.bam >null 2>&1")
            print("Calculating coverage....")
            os.system(installationDirectory+"src/conda/bin/samtools depth alignment_sorted.bam >coverage.txt 2>&1")

            covValues = {}
            print("Analyzing coverage data....")
            cov = open("coverage.txt")
            while True:
                covLine = cov.readline().rstrip()
                if not covLine:
                    break
                covFields = covLine.split("\t")
                if not covFields[0] in covValues:
                    covValues[covFields[0]] = []
                    covValues[covFields[0]].append(float(covFields[2]))
                else:
                    covValues[covFields[0]].append(float(covFields[2]))
            cov.close()
            bestCoverage = 0
            bestMatch = ""
            for item in covValues:
                print(item+"\t"+str(len(covValues[item]))+"\t"+str(nm.mean(covValues[item])))
                if nm.mean(covValues[item]) > bestCoverage:
                    bestCoverage = nm.mean(covValues[item])
                    bestMatch = item


            
            if not bestMatch == '':
                print("Genome",bestMatch,"seems to be the best match to fill the gap")
                print("The bit length is",len(foundSequences[bestMatch]),"and it is covered with an average coverage",bestCoverage,"bases")
        


                #get coverage for bestMatch
                cov = open("coverage.txt")
                coveredPositions = []
                while True:
                    covLine = cov.readline().rstrip()
                    if not covLine:
                        break
                    covFields = covLine.split("\t")
                    if covFields[0] == bestMatch and int(covFields[2])>1:
                        coveredPositions.append(int(covFields[1]))
                cov.close()
                print(coveredPositions)
                coveredSequence = ""
                for a in range(len(foundSequences[bestMatch])):
                    if a in coveredPositions:
                        coveredSequence += foundSequences[bestMatch][a]
                    else:
                        coveredSequence+="J"

                print(coveredSequence)
                #[best5hits[bestMatch][0]:best3hits[bestMatch][1]]


                newSeq = seq2fill[:gapStart-200+best5hits[bestMatch][2]]+coveredSequence+seq2fill[gapEnd+best3hits[bestMatch][3]:]
                print("The new genome sequence is now",len(newSeq))
                seq2fill = newSeq
            else:
                newSeq = seq2fill[:gapStart]
                for u in range(gapEnd-gapStart):
                    newSeq+="J"
                newSeq+=seq2fill[gapEnd:]
                seq2fill = newSeq
            

            a = 0

    
outfile = open("filledGenome.fasta","w")
outfile.write(">finalScaffold\n"+seq2fill.replace("J","N")+"\n")
outfile.close()
             
