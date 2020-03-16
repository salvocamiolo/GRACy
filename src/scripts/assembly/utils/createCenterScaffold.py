from Bio import SeqIO
import sys
import os

projectName = sys.argv[1]
installationDirectory = sys.argv[2]
insertSize = sys.argv[3]
numThreads = sys.argv[4]

os.system("cp "+installationDirectory+"data/merlinReference/hcmv_genome.fasta .")
os.system("cp "+installationDirectory+"src/scripts/assembly/utils/genome_recepie.rcp  .")
os.system("cp "+installationDirectory+"src/scripts/assembly/utils/center_recepie.rcp .")


#new bit
os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build merlinGenome_190000_200000_f.txt reference")
os.system(installationDirectory+"src/conda/bin/bowtie2  -x reference -1 ../1_cleanReads/"+projectName+"_hq_1.fastq -2 ../1_cleanReads/"+projectName+"_hq_2.fastq -S alignment.sam -p "+numThreads)
os.system(installationDirectory+"src/conda/bin/samtools view -h -Sb alignment.sam >alignment.bam")
os.system(installationDirectory+"src/conda/bin/samtools view -F 12 -b alignment.bam >botMapped.bam")
os.system(installationDirectory+"src/conda/bin/samtools view -f 4 -F 8 -b alignment.bam >oneMapped.bam")
os.system(installationDirectory+"src/conda/bin/samtools view -f 8 -F 4 -b alignment.bam >twoMapped.bam")
os.system(installationDirectory+"src/conda/bin/samtools merge all.bam botMapped.bam oneMapped.bam twoMapped.bam")
os.system(installationDirectory+"src/conda/bin/bam2fastq -o reads#.fastq all.bam")
os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/spades.py -1 reads_1.fastq -2 reads_2.fastq  --cov-cutoff auto --careful -k 51,61,71 -o centerScaffold")
#Following command replaced by the following
os.system(installationDirectory+"src/conda2/bin/python "+installationDirectory+"src/scripts/assembly/utils/scaffold_builder.py "+installationDirectory+" -q ./centerScaffold/scaffolds.fasta -r merlinGenome_190000_200000_f.txt -p sb2 ")
#os.system(installationDirectory+"resources/Ragout/bin/ragout --overwrite center_recepie.rcp")




#os.system("scaffold_builder_v2.py -q ../2_spadesAssembly/scaffolds.fasta -r merlinGenome_190000_200000_f.txt -p sb2 >null 2>&1")
#Check the sb2_scaffold file for anomalous bases
outfile = open("temp.fasta","w")
for seq_record in SeqIO.parse("sb2_Scaffold.fasta","fasta"):
    newSeq = ""
    sequence = str(seq_record.seq)
    for a in range(len(sequence)):
        if not sequence[a] =="A" and not sequence[a] =="T" and not sequence[a] =="G" and not sequence[a] =="C" and  not sequence[a] =="a" and not sequence[a] =="t" and not sequence[a] =="g" and not sequence[a] =="c":
            newSeq += "N"
        else:
            newSeq += sequence[a]
    outfile.write(">"+str(seq_record.id)+"\n"+newSeq+"\n")
outfile.close()
os.system("mv temp.fasta sb2_Scaffold.fasta")


os.system(installationDirectory+"src/conda/bin/lastz merlinGenome_190000_200000_f.txt sb2_Scaffold.fasta --format=general:name1,start1,end1,name2,start2,end2,identity,score --ydrop=50000 >lastzOutput.txt")

bestAlignment = ["scaffold",1,1,1,1,0]

infile = open("lastzOutput.txt")
infile.readline()
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if int(fields[8])>bestAlignment[5]:
        bestAlignment = [fields[3] , int(fields[1]),int(fields[2]), int(fields[4]),int(fields[5]),int(fields[8])]


print(bestAlignment)
print("The center scaffold has length",str(bestAlignment[2]-bestAlignment[1]))
#print "The best alignment is ",bestAlignment
#Add the center repetitive region to the contigs files
os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py sb2_Scaffold.fasta "+bestAlignment[0]+" "+str(bestAlignment[3])+" "+str(bestAlignment[4])+" f")
#os.system("cat ../2_spadesAssembly/scaffolds.fasta "+bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt >scaffolds.fasta ")

#Next command is changed with the following
#os.system(installationDirectory+"resources/scaffold_builder_v2.py -q ../2_spadesAssembly/scaffolds.fasta -r /home3/scc20x/hcmvReference/hcmv_genome.fasta -p sb >null 2>&1")
os.system(installationDirectory+"src/conda2/bin/ragout --overwrite genome_recepie.rcp")


for seq_record in SeqIO.parse("./ragout-out/scaffolds_scaffolds.fasta","fasta"):
    mainScaffold = str(seq_record.seq)
    break

print("The main scaffold has length",len(mainScaffold))



#Get the central repetitive region into the variable centralScaffold
for seq_record in SeqIO.parse(bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt","fasta"):
    centerScaffold = str(seq_record.seq)

#for seq_record in SeqIO.parse("sb2_Scaffold.fasta","fasta"):
#    print str(seq_record.id)
#    if str(seq_record.id) == "Scaffold_1":
#        centerScaffold = (str(seq_record.seq))[:10000]

print(len(centerScaffold))

#Check for the presence of N stretches and, if present, run gapfiller:
if "N" in centerScaffold:
    print("Ns present in the center sequence. Now running gapfiller....")
    gfFile = open("gapfillerlib.txt","w")
    gfFile.write("lib1 bwa ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq "+ insertSize+" 0.25 FR")
    gfFile.close()
    centerScaffoldFile = bestAlignment[0]+"_"+str(bestAlignment[3])+"_"+str(bestAlignment[4])+"_f.txt"
    os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/GapFiller -q "+installationDirectory+" -l gapfillerlib.txt -s "+centerScaffoldFile+" -T "+numThreads)
    os.system("cp ./standard_output/standard_output.gapfilled.final.fa "+centerScaffoldFile)




#Get the 5' and 3' recombination coordinates for the central repetitive region centralScaffold
recombinationStart = 0
recombinationEnd = 0
foundPoint = False
for a in range(500):
    if not mainScaffold[10000:].find(centerScaffold[a:a+100]) == -1:
        recombinationStart = mainScaffold[10000:].find(centerScaffold[a:a+100])
        start2 = a
        foundPoint = True

for a in range(500):
    if not mainScaffold[10000:].find(centerScaffold[-100-a:-a]) == -1:
        recombinationEnd = mainScaffold[10000:].find(centerScaffold[-100-a:-a])
        end2=a
        foundPoint = True

print("Recombination start and end for central repetitive region are ",recombinationStart+10000,recombinationEnd+10000)
if not recombinationStart==0 and not recombinationEnd==0:
    newScaffold = mainScaffold[:recombinationStart+10000]+centerScaffold[start2:len(centerScaffold)-100-end2]+mainScaffold[recombinationEnd+10000:]
else:
    newScaffold = mainScaffold


outfile = open("newFinalScaffold.fasta","w")
outfile.write(">finalScaffold\n"+newScaffold+"\n")
outfile.close()
