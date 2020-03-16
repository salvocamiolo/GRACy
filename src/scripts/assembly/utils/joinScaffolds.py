#This software elong the sequences in a fasta file by using the alignment
#of reads from the same organism and exploting the unmapped reads whose
#mate is mapped on the sequence.
#Usage:
#python elong.py projectNAme pathToRead1.fastq pathToRead2.fastq sequencesToElong.fasta numCycles


#Version 0.1



import os
import sys
from Bio import SeqIO
import biomodule

# Import pairwise2 module
from Bio import pairwise2
# Import format_alignment method
from Bio.pairwise2 import format_alignment


projectName =sys.argv[1]
read1 = sys.argv[2]
read2 = sys.argv[3]
sequenceToElong = sys.argv[4]
sequenceToElongOrientation = sys.argv[5]
sequenceToReach = sys.argv[6]
sequenceToReachOrientation = sys.argv[7]

installationDirectory = sys.argv[8]
numThreads = sys.argv[9]

def fuseSequences2(s1,s2):

	start = open("s1.fasta","w")
	start.write(">s1\n"+s1+"\n")
	start.close()
	toFuse = open("s2.fasta","w")
	toFuse.write(">s2\n"+s2+"\n")
	toFuse.close()
	os.system(installationDirectory+"src/conda/bin/makeblastdb -dbtype nucl -in s1.fasta >null 2>&1")
	os.system(installationDirectory+"src/conda/bin/blastn -query s2.fasta -db s1.fasta -outfmt 6 -dust no -soft_masking false -task blastn -out outputBlast.txt >null 2>&1")
	blastFile = open("outputBlast.txt")

	downstreamAlignment = []
	lastNucl = 0

	while True:
		line = blastFile.readline().rstrip()
		if not line:
			break
		fields = line.split("\t")
		if int(fields[9]) > lastNucl and ( int(fields[9])-int(fields[8]) )>30 and float(fields[10]) <0.0001:
			downstreamAlignment = fields
			lastNucl = int(fields[9])

	if len(downstreamAlignment) >0:
		newSequence = s1[:int(downstreamAlignment[9])]+s2[int(downstreamAlignment[7]):]
		blastFile.close()
		return newSequence
	else:
		blastFile.close()
		return ""


def greedyElongation(seq):
	unmappedSequences = {}
	reference = str(seq)
	#os.system("mkdir tempor")
	for seq_record in SeqIO.parse("unmapped.fasta","fasta"):
			locus = str(seq_record.id)
			if not locus in unmappedSequences:
				unmappedSequences[locus] = str(seq_record.seq)
	numElong = 0
	while True:
		#print "Perform the blast of the unmapped sequences...."
		toAssemble = open("toAssemble.fasta","w")
		toAssemble.write(">toElong\n"+reference[-200:]+"\n")
		os.system(installationDirectory+"src/conda/bin/makeblastdb -in toElong.fasta -dbtype nucl >null 2>&1")
		os.system(installationDirectory+"src/conda/bin/blastn -query unmapped.fasta -db toElong.fasta -outfmt 6 -num_threads 8 -dust no -soft_masking false -out outputBlast.txt  >null 2>&1")
		#print "Done"

		#print "Fill the toAssemble file"
		blastFile = open("outputBlast.txt")
		while True:
			line = blastFile.readline().rstrip()
			if not line:
				break
			fields = line.split("\t")
			if (int(fields[8]) > (len(reference)-100) or int(fields[9]) > (len(reference)-100)) and abs(int(fields[9])-int(fields[8])) > 40 and float(fields[2])>97 and fields[5] == "0":

				if int(fields[8]) > int(fields[9]):
					toAssemble.write(">"+fields[0]+"\n"+biomodule.reverseComplement(unmappedSequences[fields[0]])+"\n")
				else:
					toAssemble.write(">"+fields[0]+"\n"+unmappedSequences[fields[0]]+"\n")

		toAssemble.close()

		print("Perform second phrap assembly step.......")

		os.system(installationDirectory+"src/conda/bin/cap3 toAssemble.fasta > cap3Assembly 2>null")

		numElong += 1

		longestScaffold = ""
		for seq_record in SeqIO.parse("toAssemble.fasta.cap.contigs","fasta"):
			if len(str(seq_record.seq)) >= len(longestScaffold):
				longestScaffold = str(seq_record.seq)

		 #Check whether the produced elonged scaffold is on the right orientation
		tempScaffold = open("tempScaffold.fasta","w")
		tempScaffold.write(">tempScaffold\n"+longestScaffold+"\n")
		tempScaffold.close()
		tempScaffold = open("tempScaffold_query.fasta","w")
		tempScaffold.write(">reference\n"+reference[-100:]+"\n")
		tempScaffold.close()

		os.system(installationDirectory+"src/conda/bin/makeblastdb -in tempScaffold.fasta -dbtype nucl >null 2>&1")
		os.system(installationDirectory+"src/conda/bin/blastn -query tempScaffold_query.fasta -db tempScaffold.fasta -outfmt 6 -num_threads 10 -dust no -soft_masking false -out tempScaffold_outputBlast.txt >null 2>&1")
		tempScaffold = open("tempScaffold_outputBlast.txt")
		line = tempScaffold.readline().rstrip()
		fieldBlast = line.split("\t")
		if len(fieldBlast)>2:
			if int(fieldBlast[8]) > int(fieldBlast[9]):
				longestScaffold = biomodule.reverseComplement(longestScaffold)
			tempScaffold.close()

		else:
			print("WARNING! EXTENSION STOPPED!!")
			js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
			js.write(">joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+reference)
			js.close()
			exit()


		sc = fuseSequences2(reference,longestScaffold)
		longestScaffold = sc
		print("Elonged sequence has now a size of",len(longestScaffold),"nucleotides")

		if len(longestScaffold) <= len(reference):

			return longestScaffold
		else:
			reference = longestScaffold
			toElong = open("toElong.fasta","w")
			toElong.write(">toElong\n"+reference+"\n")
			toElong.close()



sequences = {}
for seq_record in SeqIO.parse(sequenceToElong,"fasta"):
	startingSeq = str(seq_record.seq)
	id1 = str(seq_record.id)
	if sequenceToElongOrientation == "r":
		startingSeq = biomodule.reverseComplement(startingSeq)

for seq_record in SeqIO.parse(sequenceToReach,"fasta"):
	terminiSeq = str(seq_record.seq)
	id2 = str(seq_record.id)
	if sequenceToReachOrientation == "r":
		terminiSeq = biomodule.reverseComplement(terminiSeq)

termfile = open("termini.fasta","w")
termfile.write(">termini\n"+terminiSeq[:500]+"\n")
termfile.close()

toElong = open("toElong.fasta","w")
toElong.write(">toElong\n"+startingSeq[-1800:-300]+"\n")
toElong.close()


numCycle = 0
unmappedContigs = 0
while True:

	for seq_record in SeqIO.parse("toElong.fasta","fasta"):
		reference = str(seq_record.seq)

	print("Indexing.....")
	os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build  toElong.fasta toElong >null 2>&1")
	print("Aligning.....")
	os.system(installationDirectory+"src/conda/bin/bowtie2 --local --very-fast-local -k 3 -p "+numThreads+"  -x toElong -1 "+ read1 + " -2 " + read2 + " -S "+ projectName+".sam >null 2>&1" )

	os.system(installationDirectory+"src/conda/bin/samtools view -bS -F 8 -f 4 "+ projectName+".sam > unmapped.bam 2>null")
	os.system(installationDirectory+"src/conda/bin/bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq >null 2>&1")
	os.system(installationDirectory+"src/conda/bin/fastq_to_fasta -i unmapped.fastq -o unmapped.fasta -Q33 >null 2>&1")

	print("Clustering unmapped reads.....")
	os.system(installationDirectory+"src/conda/bin/cd-hit-est  -d 0  -i unmapped.fasta -o unmapped_cdhit.fasta >null 2>&1")

	cdhitfile = open("unmapped_cdhit.fasta.clstr")
	#Select only high representated unmapped reads

	unmapSeq = {}
	for seq_record in SeqIO.parse("unmapped.fasta","fasta"):
		locus = str(seq_record.id)
		if not locus in unmapSeq:
			unmapSeq[locus] = str(seq_record.seq)


	clusters ={}
	numSeqInCluster = {}

	numCluster = 0
	line = cdhitfile.readline().rstrip()
	if not line:
		print("")
	else:
		while True:
			clusterName = "cluster"+str(numCluster)
			if not clusterName in clusters:
				clusters[clusterName] = ""
				numSeqInCluster[clusterName] = 0
			line = "start"
			while not line[0] == ">":
				line = cdhitfile.readline().rstrip()
				if not line:
					break
				if "*" in line:
					clusters[clusterName] = ((line.split(">"))[1].split("..."))[0]
					numSeqInCluster["cluster"+str(numCluster)] += 1
				else:
					numSeqInCluster["cluster"+str(numCluster)] += 1
				if line[0] == ">":
					numCluster +=1
			if not line:
				break


	outcdhitfile = open("tempCdhitFile","w")
	for item in numSeqInCluster:
		if numSeqInCluster[item]>10:
			outcdhitfile.write(">"+clusters[item]+"\n"+unmapSeq[clusters[item]]+"\n")
	outcdhitfile.close()

	numSeq = 0


	if unmappedContigs == 1:
		os.system("cat unmapped.fasta unmapped.fasta.contigs2 >tmp1 ; mv tmp1 unmapped.fasta")
	print("Assembling unmapped reads......")


	os.system(installationDirectory+"src/conda/bin/cap3 unmapped.fasta > cap3Assembly 2>null")
	unmappedContigsFile = open("unmapped.fasta.contigs2","w")

	unmappedContigs = 0
	for seq_record in SeqIO.parse("unmapped.fasta.cap.contigs","fasta"):

		tempScaffold = open("tempScaffold.fasta","w")
		tempScaffold.write(">tempScaffold\n"+str(seq_record.seq)+"\n")
		tempScaffold.close()
		os.system(installationDirectory+"src/conda/bin/makeblastdb -in tempScaffold.fasta -dbtype nucl >null 2>&1")
		os.system(installationDirectory+"src/conda/bin/blastn -query toElong.fasta -db tempScaffold.fasta -dust no -soft_masking false -outfmt 6 -out tempScaffold_outputBlast.txt >null 2>&1")
		tempScaffold = open("tempScaffold_outputBlast.txt")

		line = tempScaffold.readline().rstrip()
		fieldBlast = line.split("\t")
		if len(fieldBlast) >2:
			if abs(int(fieldBlast[8]) - int(fieldBlast[9])) >300 and int(fieldBlast[6]) < (len(reference)-200):
				print("Redundant scaffold")
			else:
				unmappedContigsFile.write(">"+str(seq_record.id)+"\n"+str(seq_record.seq)+"\n")
				unmappedContigs = 1
		else:
			unmappedContigsFile.write(">"+str(seq_record.id)+"\n"+str(seq_record.seq)+"\n")
			unmappedContigs = 1



		tempScaffold.close()
	unmappedContigsFile.close()


	alignment = open(projectName+".sam")
	softClipped = open("softclippedReads","w")
	numSoftClipped = 0
	while True:
		line = alignment.readline().rstrip()
		if not line:
			break
		fields = line.split("\t")
		if len(fields)>=6 and not fields[3]=="*":
			if fields[5][-1:] == "S" and int(fields[3]) > (len(reference)-80) and fields[5].count('S') == 1 and not "I" in fields[5] and not "D" in fields[5] and not fields[7]=="*":
				splitCigar = fields[5].split("M")
				if len(splitCigar) == 2 and int(splitCigar[0])>40:# and (int(fields[3])+int(splitCigar[0])==len(reference)+1):
					splitCigar2 = splitCigar[1].split("S")
					if len(splitCigar2) == 2 and int(splitCigar2[0])>40:
						softClipped.write(">"+fields[0]+"_softClipped\n"+fields[9]+"\n")
						numSoftClipped += 1
						if numSoftClipped >100:
							break

	softClipped.close()
	os.system("cat unmapped.fasta softclippedReads>tmp1; mv tmp1 unmapped.fasta")



	longestScaffold = greedyElongation(reference)
	print(longestScaffold)
	if len(longestScaffold) <= len(reference):
		print("WARNING! EXTENSION STOPPED DUE TO NO ELONGATION!!")
		js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
		js.write(">joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+reference)
		js.close()
		exit()

	if len(longestScaffold) > 10000:
		print("WARNING! EXTENSION STOPPED DUE TO NEVER ENDING ELONGATION!!")
		js = open("joined_W_WARNING_"+sequenceToElong+"_"+sequenceToReach,"w")
		js.write(">joined_"+sequenceToElong+"_"+sequenceToReach+"_NeverEndingElongation\n"+startingSeq+"\n")
		js.close()
		exit()





	#Check whether the produced elonged scaffold is on the right orientation
	tempScaffold = open("tempScaffold.fasta","w")
	tempScaffold.write(">tempScaffold\n"+longestScaffold+"\n")
	tempScaffold.close()
	os.system(installationDirectory+"src/conda/bin/makeblastdb -in tempScaffold.fasta -dbtype nucl >null 2>&1")
	os.system(installationDirectory+"src/conda/bin/blastn -query toElong.fasta -db tempScaffold.fasta -dust no -soft_masking false -outfmt 6 -out tempScaffold_outputBlast.txt >null 2>&1")
	tempScaffold = open("tempScaffold_outputBlast.txt")
	line = tempScaffold.readline().rstrip()
	fieldBlast = line.split("\t")
	if int(fieldBlast[8]) > int(fieldBlast[9]):
		longestScaffold = biomodule.reverseComplement(longestScaffold)
	tempScaffold.close()




	#check if the elonged sequence reached the specified termini portion
	fusedTermini = fuseSequences2(longestScaffold,terminiSeq[:500])
	if not fusedTermini == "":
		js = open("joined_"+sequenceToElong+"_"+sequenceToReach,"w")
		js.write(">joined_"+sequenceToElong+"_"+sequenceToReach+"\n"+startingSeq[:-1800]+fuseSequences2(longestScaffold,terminiSeq))
		js.close()
		exit()
