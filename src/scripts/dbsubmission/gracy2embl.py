import os
import sys
import argparse
from Bio import SeqIO



parser = argparse.ArgumentParser(description="This tool converts the GRACy annotation and assembled sequence in EMBL bank format")
parser.add_argument("-gff","--annotation",required=True, help="The GFF file produced by GRACy")
parser.add_argument("-prot","--proteinFile",required=True, help="The multifasta formatted protein file produced by GRACy")
parser.add_argument("-i","--info",required=True, help="The annotation info file in txt format")
parser.add_argument("-g","--genomeSequence",required=True, help="The assembled genomic sequence")
parser.add_argument("-o","--outputFile",required=True,help="Output file name")
args = vars(parser.parse_args())
gffFile = args['annotation']
protFile = args['proteinFile']
infoFile = args['info']
genomeFile = args['genomeSequence']
outputFile = args['outputFile']


gene2report = ["RL1","RNA2.7","RL5A","RL6","RNA1.2","RL8A","RL9A","RL10","RL11","RL12","RL13","UL1","UL2","UL4","UL5","UL6","UL8","UL7","UL9","UL10","UL11","UL13","UL14","UL15A","UL16","UL17","UL18","UL19","UL20","UL21A","UL22A","UL23","UL24","UL25","UL26","UL27","UL29","UL30","UL30A","UL31","UL32","UL33","UL34","UL35","UL36","UL37","UL38","UL40","UL41A","UL42","UL43","UL44","UL45","UL46","UL47","UL48","UL48A","UL49","UL50","UL51","UL52","UL53","UL54","UL55","UL56","UL57","RNA4.9","UL69","UL70","UL71","UL72","UL73","UL74","UL74A","UL75","UL76","UL77","UL78","UL79","UL80","UL80.5","UL82","UL83","UL84","UL85","UL86","UL87","UL88","UL89","UL91","UL92","UL93","UL94","UL95","UL96","UL97","UL98","UL99","UL100","UL102","UL103","UL104","UL105","RNA5.0","UL111A","UL112","UL114","UL115","UL116","UL117","UL119","UL120","UL121","UL122","UL123","UL124","UL128","UL130","UL131A","UL132","UL148","UL147A","UL147","UL146","UL145","UL144","UL142","UL141","UL140","UL139","UL138","UL136","UL135","UL133","UL148A","UL148B","UL148C","UL148D","UL150","UL150A","IRS1","US1","US2","US3","US6","US7","US8","US9","US10","US11","US12","US13","US14","US15","US16","US17","US18","US19","US20","US21","US22","US23","US24","US26","US27","US28","US29","US30","US31","US32","US33A","US34","US34A","TRS1"]
productNames = {"RL1":"protein RL1","RL8A":"protein RL8A","RL9A":"protein RL9A","RL10":"envelope glycoprotein RL10","RL11":"membrane glycoprotein RL11","RL12":"membrane protein RL12","RL13":"membrane protein RL13","UL1":"membrane protein UL1","UL2":"protein UL2","UL4":"envelope glycoprotein UL4","UL5":"protein UL5","UL6":"membrane protein UL6","UL8":"membrane glycoprotein UL8","UL7":"membrane glycoprotein UL7","UL9":"membrane glycoprotein UL9","UL10":"membrane protein UL10","UL11":"membrane glycoprotein UL11","UL13":"protein UL13","UL14":"membrane protein UL14","UL15A":"protein UL15A","UL16":"membrane glycoprotein UL16","UL17":"protein UL17","UL18":"membrane glycoprotein UL18","UL19":"protein UL19","UL20":"membrane protein UL20","UL21A":"protein UL21A","UL22A":"glycoprotein UL22A","UL23":"tegument protein UL23","UL24":"tegument protein UL24","UL25":"tegument protein UL25","UL26":"tegument protein UL26","UL27":"protein UL27","UL29":"protein UL29","UL30":"protein UL30","UL31":"protein UL31","UL32":"tegument protein pp150","UL33":"envelope glycoprotein UL33","UL34":"protein UL34","UL35":"tegument protein UL35","UL36":"tegument protein vICA","UL37":"envelope glycoprotein UL37","UL38":"protein UL38","UL40":"membrane glycoprotein UL40","UL41A":"protein UL41A","UL42":"protein UL42","UL43":"tegument protein UL43","UL44":"DNA polymerase processivity subunit","UL45":"ribonucleotide reductase subunit 1","UL46":"capsid triplex subunit 1","UL47":"tegument protein UL37","UL48":"large tegument protein","UL48A":"small capsid protein","UL49":"protein UL49","UL50":"nuclear egress membrane protein","UL51":"DNA packaging protein UL33","UL52":"DNA packaging protein UL32","UL53":"nuclear egress lamina protein","UL54":"DNA polymerase catalytic subunit","UL55":"envelope glycoprotein B","UL56":"DNA packaging terminase subunit 2","UL57":"single-stranded DNA-binding protein","UL69":"multifunctional expression regulator","UL70":"helicase-primase primase subunit","UL71":"tegument protein UL51","UL72":"deoxyuridine triphosphatase","UL73":"envelope glycoprotein N","UL74":"envelope glycoprotein O","UL74A":"envelope glycoprotein 24","UL75":"envelope glycoprotein H","UL76":"nuclear protein UL24","UL77":"DNA packaging tegument protein UL25","UL78":"envelope protein UL78","UL79":"protein UL79","UL80":"capsid maturation protease","UL80.5":"capsid scaffold protein","UL82":"tegument protein pp71","UL83":"tegument protein pp65","UL84":"protein UL84","UL85":"capsid triplex subunit 2","UL86":"major capsid protein","UL87":"protein UL87","UL88":"tegument protein UL88","UL89":"DNA packaging terminase subunit 1","UL91":"protein UL91","UL92":"protein UL92","UL93":"DNA packaging tegument protein UL17","UL94":"tegument protein UL16","UL95":"protein UL95","UL96":"tegument protein UL14","UL97":"tegument serine/threonine protein kinase","UL98":"deoxyribonuclease","UL99":"myristylated tegument protein","UL100":"envelope glycoprotein M","UL102":"helicase-primase subunit","UL103":"tegument protein UL7","UL104":"capsid portal protein","UL105":"helicase-primase helicase subunit","UL111A":"interleukin-10","UL112":"protein UL112","UL114":"uracil-DNA glycosylase","UL115":"envelope glycoprotein L","UL116":"protein UL116","UL117":"protein UL117","UL119":"membrane glycoprotein UL119","UL120":"membrane protein UL120","UL121":"membrane protein UL121","UL122":"regulatory protein IE2","UL123":"regulatory protein IE1","UL124":"membrane protein UL124","UL128":"envelope protein UL128","UL130":"envelope glycoprotein UL130","UL131A":"envelope protein UL131A","UL132":"envelope glycoprotein UL132","UL148":"membrane protein UL148","UL147A":"membrane protein UL147A","UL147":"chemokine vCXCL2","UL146":"chemokine vCXCL1","UL145":"protein UL145","UL144":"membrane glycoprotein UL144","UL142":"membrane glycoprotein UL142","UL141":"membrane glycoprotein UL141","UL140":"protein UL140","UL139":"membrane glycoprotein UL139","UL138":"protein UL138","UL136":"protein UL136","UL135":"protein UL135","UL133":"protein UL133","UL148A":"protein UL148A","UL148B":"protein UL148B","UL148C":"protein UL148C","UL148D":"protein UL148D","UL150":"protein UL150","UL150A":"protein UL150A","IRS1":"tegument protein IRS1","US1":"protein US1","US2":"membrane glycoprotein US2","US3":"membrane glycoprotein US3","US6":"membrane glycoprotein US6","US7":"membrane glycoprotein US7","US8":"membrane glycoprotein US8","US9":"membrane glycoprotein US9","US10":"membrane glycoprotein US10","US11":"membrane glycoprotein US11","US12":"membrane protein US12","US13":"membrane protein US13","US14":"membrane protein US14","US15":"membrane protein US15","US16":"membrane protein US16","US17":"membrane protein US17","US18":"membrane protein US18","US19":"membrane protein US19","US20":"membrane protein US20","US21":"membrane protein US21","US22":"tegument protein US22","US23":"protein US23","US24":"tegument protein US24","US26":"protein US26","US27":"envelope glycoprotein US27","US28":"envelope protein US28","US29":"membrane protein US29","US30":"membrane protein US30","US31":"protein US31","US32":"protein US32","US33A":"protein US33A","US34":"protein US34","US34A":"protein US34A","TRS1":"tegument protein TRS1"}

#collect genomic sequence
for seq_record in SeqIO.parse(genomeFile,"fasta"):
    genomeSeq = str(seq_record.seq)

#Collect protein sequences
protSeqs = {}
protStrand = {}
for seq_record in SeqIO.parse(protFile,"fasta"):
    protID = str(seq_record.id)
    sequence = str(seq_record.seq)
    if not protID in protSeqs:
        protSeqs[protID] = sequence.replace("*","")
    if not protID in protStrand:
        protStrand[protID] = ((str(seq_record.description)).split(" "))[1]

#Collect the gene coordinates
infile = open(gffFile)
geneCoord = {}
while True:
    line = infile.readline().rstrip()
    if not line:
        break
    fields = line.split("\t")
    if fields[2] == "CDS":
        geneName = (fields[8].split("Product="))[-1]
        if not geneName in geneCoord:
            geneCoord[geneName] = ""
        if len(geneCoord[geneName])>0:
            geneCoord[geneName]+=","
        if geneName == "UL30A" and genomeSeq[int(fields[4]):int(fields[4])+3] == "CGT":
            geneCoord[geneName]+=fields[3]+".."+str(int(fields[4])+3)
            protSeqs[geneName] = "T"+protSeqs[geneName]
        else:
            geneCoord[geneName]+=fields[3]+".."+fields[4]

infile.close()


#check coordinates for missing stop codon
for gene in gene2report:
    if gene in geneCoord:
        geneStart = int((((geneCoord[gene].split(","))[0]).split(".."))[0])-1
        geneEnd = int((((geneCoord[gene].split(","))[-1]).split(".."))[-1])

        if protStrand[gene] == "+":
            if not genomeSeq[geneEnd-3:geneEnd] == "TAG" and not genomeSeq[geneEnd-3:geneEnd] == "TGA" and not genomeSeq[geneEnd-3:geneEnd]== "TAA":
                print(gene,genomeSeq[geneEnd-3:geneEnd],genomeSeq[geneEnd:geneEnd+3])
                if genomeSeq[geneEnd:geneEnd+3] == "TAA" or genomeSeq[geneEnd:geneEnd+3] == "TGA" or genomeSeq[geneEnd:geneEnd+3] == "TAG":
                    geneCoordList = geneCoord[gene].split(",")
                    newGeneCoordList = []
                    for item in range(len(geneCoordList)):
                        if not item == len(geneCoordList)-1:
                            newGeneCoordList.append(geneCoordList[item])
                        else:
                            newGeneCoordList.append((geneCoordList[item].split(".."))[0]+".."+str(int((geneCoordList[item].split(".."))[-1])+3))

                    geneCoord[gene]= ",".join(newGeneCoordList)


        else:
            if not genomeSeq[geneStart:geneStart+3] == "CTA" and not genomeSeq[geneStart:geneStart+3] == "TCA" and not genomeSeq[geneStart:geneStart+3] == "TTA":
                print(gene,genomeSeq[geneStart:geneStart+3],genomeSeq[geneStart-3:geneStart])
                if genomeSeq[geneStart-3:geneStart] == "CTA" or genomeSeq[geneStart-3:geneStart] == "TCA" or genomeSeq[geneStart-3:geneStart] == "TTA":
                
                    geneCoordList = geneCoord[gene].split(",")
                    newGeneCoordList = []
                    for item in range(len(geneCoordList)):
                        if not item == 0:
                            newGeneCoordList.append(geneCoordList[item])
                        else:
                            newGeneCoordList.append( str(int((geneCoordList[item].split(".."))[0])-3)+".."+(geneCoordList[item].split(".."))[-1])

                    print(geneCoord[gene])
                    print(",".join(newGeneCoordList))

                    geneCoord[gene]= ",".join(newGeneCoordList)

                


#Collecting data from info file
infile = open(infoFile)
generalInfo  = ((infile.readline().rstrip()).split("\t"))[1]
assemblyMethod = ((infile.readline().rstrip()).split("\t"))[1]
coverage = ((infile.readline().rstrip()).split("\t"))[1]
sequencingTechnology = ((infile.readline().rstrip()).split("\t"))[1]
locusTag = ((infile.readline().rstrip()).split("\t"))[1]
strainName = ((infile.readline().rstrip()).split("\t"))[1]
isolationSource = ((infile.readline().rstrip()).split("\t"))[1]
country = ((infile.readline().rstrip()).split("\t"))[1]
collectionDate = ((infile.readline().rstrip()).split("\t"))[1]
notes = ((infile.readline().rstrip()).split("\t"))[1]
infile.close()

#writing the EMBL fortmat
outfile = open(outputFile,"w")
outfile.write("ID   "+strainName+"; SV 0; linear; genomic DNA; STD; VRL; "+str(len(genomeSeq))+" BP.\n")
outfile.write("XX\nAC   ;\nXX\n")
outfile.write("DE   "+generalInfo+"\nXX\nKW   .\nXX\nOS   Human herpesvirus 5 (Human cytomegalovirus)\nOC   Viruses; dsDNA viruses, no RNA stage; Herpesviridae; Betaherpesvirinae;\nOC   Cytomegalovirus.\nXX\nCC   ##Assembly-Data-START##\n")
outfile.write("CC   Assembly Method       :: "+assemblyMethod+"\n")
outfile.write("CC   Coverage              :: "+coverage+"\n")
outfile.write("CC   Sequencing Technology :: "+sequencingTechnology+"\n")
outfile.write("CC   ##Assembly-Data-END##\n")
outfile.write("AC * _"+locusTag+"\nXX\n")
outfile.write("FH   Key             Location/Qualifiers\nFH\n")
outfile.write("FT   source          1.."+str(len(genomeSeq))+"\n")
outfile.write("FT                   /organism=\"Human herpesvirus 5\""+"/mol_type=\"genomic DNA\"\n")
outfile.write("FT                   /strain=\""+strainName+"\"\n")
outfile.write("FT                   /segment=\"whole_genome\"\n"+"FT                   /isolation_source=\""+isolationSource+"\"\n")
outfile.write("FT                   /host=\"Homo sapiens\"\n"+"FT                   /country=\""+country+"\"\n")
outfile.write("FT                   /collection_date=\""+collectionDate+"\"\n")
outfile.write("FT                   /note=\""+notes+"\"\n")
outfile.write("FT                   acronym: HCMV; acronym: \"HHV-5\"\n")

for gene in gene2report:
    if gene in protSeqs:
        print(protSeqs[gene])
        print(geneCoord[gene])
        print(protStrand[gene])
        geneStartCoord = (((geneCoord[gene].split(","))[0]).split(".."))[0]
        geneEndCoord = (((geneCoord[gene].split(","))[-1]).split(".."))[-1]
        print(geneStartCoord,geneEndCoord)
        numExons = geneCoord[gene].count("..")
        if numExons>1:
            geneCoord[gene] = "join("+geneCoord[gene]+")"
        if protStrand[gene]=="-":
            geneCoord[gene] = "complement("+geneCoord[gene]+")"
        
        if protStrand[gene] == "-":
            outfile.write("FT   gene            complement("+geneStartCoord+".."+geneEndCoord+")\n")
        else:
            outfile.write("FT   gene            "+geneStartCoord+".."+geneEndCoord+"\n")
        outfile.write("FT                   /gene=\""+gene+"\"\n")
        outfile.write("FT   CDS             "+geneCoord[gene]+"\n")
        outfile.write("FT                   /gene=\""+gene+"\"\n")
        outfile.write("FT                   /codon_start=1\n")
        ##handle not for UL30A
        if gene == "UL30A":
            outfile.write("FT                   /note=\"alternative start codon\"\n")
        if gene == "UL4" or gene == "UL6":
            outfile.write("FT                   /note=\"Multiple starting codons\"\n")
        if gene == "RL5A" or gene == "RL6":
            outfile.write("FT                   /note=\"Possible pseudogene\"\n")

        
        if gene in productNames:
            outfile.write("FT                   /product=\""+productNames[gene]+"\"\n")
        
        translationString = "FT                   /translation=\""
        if len(protSeqs[gene])<=43:
            translationString+=protSeqs[gene]
            translationString+="\n"

        else:
            translationString+=protSeqs[gene][:43]
            translationString+="\n"
            protPos = 43
            while (protPos + 57) < len(protSeqs[gene]):
                translationString+="FT                   "
                translationString+=protSeqs[gene][protPos:protPos+57]
                translationString+="\n"
                protPos = protPos+57
            translationString+="FT                   "
            translationString+=protSeqs[gene][protPos:]
            translationString+="\n"
            print(translationString)
        outfile.write(translationString)

outfile.write("SQ   Sequence "+str(len(genomeSeq))+";\n")


#write Genome sequence
seqPos = 0
while (seqPos + 60) < len(genomeSeq):
    subSeq = genomeSeq[seqPos:seqPos+60]
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]+" "+subSeq[20:30]+" "+subSeq[30:40]+" "+subSeq[40:50]+" "+subSeq[50:60]
    subSeqSpaced = "     "+subSeqSpaced
    for a in range(10-len(str(seqPos+60))):
        subSeqSpaced+=" "
    subSeqSpaced+=str(seqPos+60)
    print(subSeqSpaced)
    seqPos+=60
    outfile.write(subSeqSpaced+"\n")
subSeq = genomeSeq[seqPos:]
if len(subSeqSpaced)> 50:
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]+" "+subSeq[20:30]+" "+subSeq[30:40]+" "+subSeq[40:50]+" "+subSeq[50:60]
if len(subSeqSpaced)<= 50 and len(subSeqSpaced)>40:
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]+" "+subSeq[20:30]+" "+subSeq[30:40]+" "+subSeq[40:50]
if len(subSeqSpaced)<= 40 and len(subSeqSpaced)>30:
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]+" "+subSeq[20:30]+" "+subSeq[30:40]
if len(subSeqSpaced)<= 30 and len(subSeqSpaced)>20:
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]+" "+subSeq[20:30]
if len(subSeqSpaced)<= 20 and len(subSeqSpaced)>10:
    subSeqSpaced = subSeq[0:10]+" "+subSeq[10:20]
if len(subSeqSpaced)<= 10:
    subSeqSpaced = subSeq[0:10]
subSeqSpaced = "     "+subSeqSpaced
print(len(subSeqSpaced),len(str(genomeSeq)))
for a in range(80-len(subSeqSpaced)-len(str(len(genomeSeq)))):
        subSeqSpaced+=" "
subSeqSpaced+=str(len(genomeSeq))

outfile.write(subSeqSpaced+"\n")




outfile.write("//\n")
outfile.close()



    