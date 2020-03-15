import sys
import biomodule
from Bio import SeqIO
import os

reads1 = sys.argv[1]
reads2 = sys.argv[2]

start = sys.argv[3]
start_o = sys.argv[4]
end = sys.argv[5]
end_o = sys.argv[6]

numCycles = int(sys.argv[7])

installationDirectory = sys.argv[8]




sequences = {}
#Creating fasta concatenated file
outfile = open("reads.fasta","w")
print("Converting reads_1.fastq into fasta")
for seq_record in SeqIO.parse(reads1,"fastq"):

    if not str(seq_record.id)+"_1" in sequences:
        sequences[str(seq_record.id)+"_1"] = str(seq_record.seq)
    outfile.write(">"+str(seq_record.id)+"_1\n"+str(seq_record.seq)+"\n")



print("Converting reads_2.fastq into fasta")
for seq_record in SeqIO.parse(reads2,"fastq"):
    if not str(seq_record.id)+"_2" in sequences:
        sequences[str(seq_record.id)+"_2"] = str(seq_record.seq)
    outfile.write(">"+str(seq_record.id)+"_2\n"+str(seq_record.seq)+"\n")
outfile.close()








def fuseSequences2(s1,s2):
    start = open("s1.fasta","w")
    start.write(">s1\n"+s1+"\n")
    start.close()
    toFuse = open("s2.fasta","w")
    toFuse.write(">s2\n"+s2+"\n")
    toFuse.close()
    os.system(installationDirectory+"src/conda/bin/makeblastdb -dbtype nucl -in s1.fasta >null 2>&1")
    os.system(installationDirectory+"src/conda/bin/blastn -query s2.fasta -db s1.fasta -outfmt 6 -task blastn  -dust no -soft_masking false -out outputBlast.txt >null 2>&1")
    blastFile = open("outputBlast.txt")

    downstreamAlignment = []
    lastNucl = 0

    while True:
        line = blastFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if int(fields[9]) > (len(s1)-35) and ( int(fields[9])-int(fields[8]) )>30 and float(fields[10]) <0.0001:
            downstreamAlignment = fields
            lastNucl = int(fields[9])

    if len(downstreamAlignment) >0:
        newSequence = s1[:int(downstreamAlignment[9])]+s2[int(downstreamAlignment[7]):]
        blastFile.close()
        return newSequence
    else:
        return ""




for seq_record in SeqIO.parse(start,"fasta"):
    startSeq = str(seq_record.seq)
    if start_o == "r":
        startSeq = biomodule.reverseComplement(startSeq)

for seq_record in SeqIO.parse(end,"fasta"):
    terminiSeq = str(seq_record.seq)
    if end_o == "r":
        terminiSeq = biomodule.reverseComplement(terminiSeq)


elongedSequence = startSeq[-700:-200]
outputSeq = open("joinScaffold_trivialSeq.fasta","w")
numCycle = 0
while True:
    bestElongation = 0
    numCycle += 1
    if numCycle == numCycles:
        outputSeq.write(">trivialSeq\n"+startSeq[:-700]+elongedSequence+"\n")
        outputSeq.close()
        exit()
    print("Performing blat cycle ",numCycle)

    for overlap in range(100,20,-20):
        print("Trying overlap",overlap,"....")
        tempFile = open("tempFasta.fasta","w")
        tempFile.write(">tempFasta\n"+elongedSequence[-overlap:]+"\n")
        tempFile.close()
        os.system(installationDirectory+"src/conda/bin/blat reads.fasta tempFasta.fasta -t=dna -q=dna -out=blast8 output.psl -fastMap >null 2>&1")

        infile = open("output.psl")
        for homology in range(100,95,-1):
            for score in range(overlap,overlap-5,-1):
                for alignmentLen in range(80,10,-20):
                    infile.seek(0)
                    foundSequences = []
                    while True:
                        line = infile.readline().rstrip()
                        if not line:
                            break
                        fields = line.split("\t")
                        if float(fields[2]) >= float(homology) and int(fields[3]) >= score:
                            readLen = len(sequences[fields[1]])
                            if int(fields[8]) < int(fields[9]): #Alignment forward
                                if readLen - int(fields[9]) >=alignmentLen:
                                    foundSequences.append(sequences[fields[1]][int(fields[9])-1:])
                            else:
                                if int(fields[9]) >=alignmentLen:
                                    foundSequences.append(biomodule.reverseComplement(sequences[fields[1]][:int(fields[9])]))


                    if len(foundSequences)>=3:
                        break
                if len(foundSequences)>=3:
                    break
            if len(foundSequences)>=3:
                break
        if len(foundSequences)>=3:
                break



    consensus = ""
    for a in range(alignmentLen):
        numA = 0
        numT = 0
        numC = 0
        numG = 0
        maxBase = 0
        calledBase = ""
        for item in foundSequences:
            if item[a] == "A":
                numA += 1
                if numA > maxBase:
                    maxBase = numA
                    calledBase = "A"
            if item[a] == "T":
                numT += 1
                if numT > maxBase:
                    maxBase = numT
                    calledBase = "T"
            if item[a] == "G":
                numG += 1
                if numG > maxBase:
                    maxBase = numG
                    calledBase = "G"
            if item[a] == "C":
                numC += 1
                if numC > maxBase:
                    maxBase = numC
                    calledBase = "C"
        consensus += calledBase


    for item in foundSequences:
        print(item)

    print("Consensus")
    print(consensus)
    if len(consensus) <2:
        outputSeq.write(">trivialSeq\n"+startSeq[:-700]+elongedSequence+"\n")
        outputSeq.close()
        exit()
    elongedSequence = elongedSequence[:-1] + consensus
    print(elongedSequence)
    print("Score",score)
    print("Homology",homology)
    print("Alignment length",alignmentLen)
    print("Number of found sequences",len(foundSequences))
