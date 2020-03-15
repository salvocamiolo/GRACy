import sys
import os
from Bio import SeqIO
import biomodule as bm


def fuseSequences2(s1, s2):
    start = open("s1.fasta", "w")
    start.write(">s1\n"+s1+"\n")
    start.close()
    toFuse = open("s2.fasta", "w")
    toFuse.write(">s2\n"+s2+"\n")
    toFuse.close()
    os.system(installationDirectory +"/resources/makeblastdb -dbtype nucl -in s1.fasta  >null 2>&1")
    os.system(installationDirectory+"/resources/blastn -query s2.fasta -db s1.fasta -outfmt 6 -dust no -soft_masking false -task blastn -out outputBlast.txt  >null 2>&1")
    blastFile = open("outputBlast.txt")

    downstreamAlignment = []
    lastNucl = 0

    while True:
        line = blastFile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        if int(fields[9]) > lastNucl and float(fields[2]) > 90.0 and (int(fields[9])-int(fields[8])) > 30 and float(fields[10]) < 0.0001:
            downstreamAlignment = fields
            lastNucl = int(fields[9])

    if len(downstreamAlignment) > 0:
        newSequence = s1[:int(downstreamAlignment[9])] + \
            s2[int(downstreamAlignment[7]):]
        blastFile.close()
        return newSequence
    else:
	blastFile.close()
        return ""

genomeToComplete = sys.argv[2]
installationDirectory =  sys.argv[1]

#attempt reconstruction 5' end
os.system(installationDirectory+"resources/lastz "+genomeToComplete+" "+installationDirectory+"/assembly/scripts/sequence_trl.fasta --format=general > outputLastz")

#get coordinates of best lastz alignment
infile = open("outputLastz")
bestScore = 0
bestPos1 = 0
bestPos2 = 0

infile.readline()
while True:
	line = infile.readline().rstrip()
	if not line:
		break
	fields = line.split("\t")
	if int(fields[0]) > bestScore:
		bestScore = int(fields[0])
		bestPos1 = fields[4]
		bestPos2 = fields[5]


#cutting the best alignment sequence from the genome to complete
os.system(installationDirectory+"/resources/extractSeqByRange.py "+genomeToComplete+" finalScaffold "+str(bestPos1)+" "+str(bestPos2)+" r")
os.system("cp "+installationDirectory+"/assembly/scripts/joinScaffolds_careful.py .")
os.system("cp "+installationDirectory+"/assembly/scripts/biomodule.py .")
os.system(installationDirectory+"/resources/extractSeqByRange.py "+genomeToComplete+" finalScaffold 1 2000 f")


os.system("python joinScaffolds_careful.py join ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq finalScaffold_1_2000_f.txt r finalScaffold_"+str(bestPos1)+"_"+str(bestPos2)+"_r.txt r  ~/Software/mySoftware/GRACy/ 8")
for seq_record in SeqIO.parse(genomeToComplete, "fasta"):
	genomeToCompleteSeq = str(seq_record.seq)
if os.path.isfile("joined_finalScaffold_1_2000_f.txt_finalScaffold_"+bestPos1+"_"+bestPos2+"_r.txt") == True:
	print "5' end successfully reconstructed!"
	
	for seq_record in SeqIO.parse("joined_finalScaffold_1_2000_f.txt_finalScaffold_"+bestPos1+"_"+bestPos2+"_r.txt", "fasta"):
		firstPortion = str(seq_record.seq)
		firstPortion = bm.reverseComplement(firstPortion)

	firtPortionReconstructed = fuseSequences2(firstPortion, genomeToCompleteSeq)
	if len(firtPortionReconstructed)>10:
		print "firstPortion successuffly joined!"
		genomeToCompleteSeq = firtPortionReconstructed
		outfile = open("newGenome1.fasta","w")
		outfile.write(">finalScaffold\n"+firtPortionReconstructed+"\n")
		outfile.close()

	else:
		print "firstPortion not joined"
		outfile = open("newGenome1.fasta", "w")
		outfile.write(">finalScaffold\n"+genomeToCompleteSeq+"\n")
		outfile.close()
	
else:
	print "5' end not reconstructed"
	outfile = open("newGenome1.fasta", "w")
	outfile.write(">finalScaffold\n"+genomeToCompleteSeq+"\n")
	outfile.close()


#attempt reconstruction of 3'end
os.system(installationDirectory+"resources/lastz newGenome1.fasta "+installationDirectory+"/assembly/scripts/sequence_irs.fasta --format=general > outputLastz")

#get coordinates of best lastz alignment
infile = open("outputLastz")
bestScore = 0
bestPos1 = 0
bestPos2 = 0

infile.readline()
while True:
	line = infile.readline().rstrip()
	if not line:
		break
	fields = line.split("\t")
	if int(fields[0]) > bestScore:
		bestScore = int(fields[0])
		bestPos1 = fields[4]
		bestPos2 = fields[5]


#cutting the best alignment sequence from the genome to complete
os.system(installationDirectory+"/resources/extractSeqByRange.py newGenome1.fasta finalScaffold "+bestPos1+" "+bestPos2+" r")
os.system("cp "+installationDirectory+"/assembly/scripts/joinScaffolds_careful.py .")
os.system("cp "+installationDirectory+"/assembly/scripts/biomodule.py .")
os.system(installationDirectory+"/resources/extractSeqByRange.py newGenome1.fasta finalScaffold "+str(len(genomeToCompleteSeq)-2000) +" "+ str(len(genomeToCompleteSeq)) + " f")

for seq_record in SeqIO.parse("newGenome1.fasta", "fasta"):
	genomeToCompleteSeq = str(seq_record.seq)
os.system("python joinScaffolds_careful.py join ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq finalScaffold_"+str(len(genomeToCompleteSeq)-2000) + "_" + str(len(genomeToCompleteSeq))+"_f.txt f finalScaffold_"+bestPos1+"_"+bestPos2+"_r.txt f  ~/Software/mySoftware/GRACy/ 8")
if os.path.isfile("joined_finalScaffold_"+str(len(genomeToCompleteSeq)-2000) + "_" + str(len(genomeToCompleteSeq))+"_f.txt_finalScaffold_"+bestPos1+"_"+bestPos2+"_r.txt") == True:
	print "3' end successfully reconstructed!"
	
	for seq_record in SeqIO.parse("joined_finalScaffold_"+str(len(genomeToCompleteSeq)-2000) + "_" + str(len(genomeToCompleteSeq))+"_f.txt_finalScaffold_"+bestPos1+"_"+bestPos2+"_r.txt","fasta"):
		lastPortion = str(seq_record.seq)

	lastPortionReconstructed = fuseSequences2(genomeToCompleteSeq, lastPortion)
	if len(lastPortionReconstructed) > 10:
		print "last portion successuffly joined!"
		genomeToCompleteSeq = lastPortionReconstructed
		outfile = open("newGenome2.fasta", "w")
		outfile.write(">finalScaffold\n"+genomeToCompleteSeq+"\n")
		outfile.close()

	else:
		print "firstPortion not joined"
		outfile = open("newGenome2.fasta", "w")
		outfile.write(">finalScaffold\n"+genomeToCompleteSeq+"\n")
		outfile.close()
	
else:
	print "5' end not reconstructed"
	outfile = open("newGenome2.fasta", "w")
	outfile.write(">finalScaffold\n"+genomeToCompleteSeq+"\n")
	outfile.close()
