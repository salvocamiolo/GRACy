import tkinter as tk
from tkinter import ttk

import tkinter
from tkinter import filedialog as tkFileDialog
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

from PIL import ImageTk, Image

import time
import os
from os import listdir
from os.path import isfile, join
import time

installationDirectory = sys.argv[1]


def vp_start_gui():
	'''Starting point when module is the main routine.'''
	global val, w, root

	root = tk.Tk()



	top = Toplevel1 (root)
	root.mainloop()



def create_Toplevel1(root, *args, **kwargs):

	'''Starting point when module is imported by another program.'''
	global w, w_win, rt
	rt = root
	w = tk.Toplevel (root)
	VirHosFilt_support.set_Tk_var()
	top = Toplevel1 (w)
	VirHosFilt_support.init(w, top, *args, **kwargs)
	return (w, top)

def destroy_Toplevel1():
	global w
	w.destroy()
	w = None

class Toplevel1:
	def __init__(self, top=None):

		def openInputFile():
			inputFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select an input file")
			self.inputFileEntry.delete(0,tk.END)
			self.inputFileEntry.insert(0,inputFile)

		def openOutputFolder():
			outputFolder = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
			self.outputFolderEntry.delete(0,tk.END)
			self.outputFolderEntry.insert(0,outputFolder)

		def exitProgram():
			exit()
		def reverseComplement(sequence):
			outSequence = ""
			accepted = ["A","C","T","G","a","c","t","g","N","n"]
			for a in range(len(sequence)-1,-1,-1):
				if sequence[a]=="A" or sequence[a]=="a":
					outSequence+="T"
				elif sequence[a]=="T" or sequence[a]=="t":
					outSequence+="A"
				elif sequence[a]=="G" or sequence[a]=="g":
					outSequence+="C"
				elif sequence[a]=="C" or sequence[a]=="c":
					outSequence+="G"
				else:
					outSequence+="N"
				

			return outSequence

		def checkCDS(cds):
			for a in range(0,len(cds)-3,+3):
				if cds[a:a+3] == "TAA" or cds[a:a+3] == "TGA" or cds[a:a+3] == "TAG":
					return "Stop codon in CDS"
			if not float(len(cds))%3.0==0:
				return "Not multiple of 3 coding sequence length"

		def mainAnnotationAlgorithm():
			if self.inputFileEntry.get() == "Please select a genomes list...." or self.outputFolderEntry.get() == "Please select an output folder....":
				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Please select a valid input file and output folder....\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()
			else:
				#**************************************************************************
				#***************** Main Annotation Algorithm Start ************************
				#**************************************************************************
				infile = open(self.inputFileEntry.get())
				outputFolder = self.outputFolderEntry.get()
				while True:
					file2Annotate = infile.readline().rstrip()
					if not file2Annotate:
						break
					file2AnnotateName = (file2Annotate.split("/"))[-1]
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Starting annotation on sample "+file2AnnotateName+"....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					genomeName = file2Annotate
					prot2map = []
					suffixName = (genomeName.split("/"))[-1]

					gffFile = open(suffixName+"_annotation.gff","w") #To change from the command line
					warnFile = open(suffixName+"_annotationWarnings.txt","w") #To change from the command line
					cdsFile = open(suffixName+"_cds.fasta","w") #To change from the command line
					protFile = open(suffixName+"_proteins.fasta","w") #To change from the command line


					os.system("rm -f "+installationDirectory+"src/scripts/annotation/proteinDB/._*")

					for seq_record in SeqIO.parse(genomeName,"fasta"):
						genomeSeq = str(seq_record.seq)
						assemblyName = str(seq_record.id)


					onlyfiles = [f for f in listdir(installationDirectory+"src/scripts/annotation/proteinDB/") if isfile(join(installationDirectory+"src/scripts/annotation/proteinDB/", f))]
					for f in onlyfiles:
						if f[-13:] == "_models.fasta":
							prot2map.append(f)

					for f in prot2map:
						numCodonRefines = 0
						#Find the best model for the gene *********************************************************
						locus = f.replace("_models.fasta","")
						protSeqs = SeqIO.to_dict(SeqIO.parse(installationDirectory+"src/scripts/annotation/proteinDB/"+f,"fasta"))
						#print "Choosing best match for protein",f
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "\nSearching best match of "+f+" for "+file2AnnotateName+"....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						os.system(installationDirectory+"src/conda/bin/makeblastdb -dbtype nucl -in "+genomeName)
						os.system(installationDirectory+"src/conda/bin/tblastn -query "+installationDirectory +
								  "src/scripts/annotation/proteinDB/"+f+" -db "+genomeName+" -outfmt 6 -max_intron_length 350000 | sort -k 12rn,12rn  > preOutputBlast.txt")
						time.sleep(1)
						os.system("head -20 preOutputBlast.txt > preOutputBlast2.txt")
						time.sleep(1)
						os.system("sort -k 3rn,3rn -k12rn,12rn preOutputBlast2.txt > outputBlast.txt")
						os.system("rm preOutputBlast.txt preOutputBlast2.txt -f ")
						blastFile = open("outputBlast.txt")
						bestCoverage = 0
						blastLine = blastFile.readline().rstrip()
						blastField = blastLine.split("\t")
						alignmentLength = int(blastField[7]) - int(blastField[6])
						bestCoverage = alignmentLength
						bestProt = blastField[0]

						for hit in range(5):
							if float(blastField[2]) > 99.8:
								break
							blastLine = blastFile.readline().rstrip()
							if not blastLine:
								break
							blastField = blastLine.split("\t")
							alignmentLength = int(blastField[7]) - int(blastField[6])
							if alignmentLength > bestCoverage:
								bestCoverage = alignmentLength
								bestProt = blastField[0]


						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Found "+bestProt+"\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						#print "Best match in the following protein:",bestProt

						sequence = str(protSeqs[bestProt].seq)
						tempFasta = open("tempFasta.fasta","w")
						tempFasta.write(">"+locus+"\n"+str(protSeqs[bestProt].seq)+"\n")
						tempFasta.close()

						#Run Exonerate on the best model *********************************************************
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Retrieving model....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --forcegtag --minintron 35 --maxintron 10000  >outputExonerate")



						#Check Exonerate output *****************************************************************
						#Check the proteins gave a match in the target genome
						exResult = open("outputExonerate")
						line = exResult.readline().rstrip()
						while not "Query range:" in line:
							line = exResult.readline().rstrip()
							if line is None:
								print("WARNING no protein found for ", locus)
								warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
						exResult.close()


						#Reconstruct Exons  **********************************************************************
						exResult = open("outputExonerate")
						while not line == "# --- START OF GFF DUMP ---":
							line = exResult.readline().rstrip()
							if line is None:
								#print "WARNING no protein found for ", locus
								warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
						for a in range(10):
							line = exResult.readline().rstrip()

						gene = {}
						exon = {}
						warnings = []
						#Collect the exonerate output
						while not line == "# --- END OF GFF DUMP ---":
							line = exResult.readline().rstrip()
							if line is None:
								#print "WARNING no protein found for ", locus
								warnings.append("Missing gene: locus "+locus+" did not provide any alignment")
							fields = line.split("\t")
							if not line == "# --- END OF GFF DUMP ---":
								if fields[2]=="gene":
									if not locus in gene:
										gene[locus] = (fields[3],fields[4],fields[6])
									else:
										if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
											gene[locus] = (fields[3],fields[4],fields[6])

								if fields[2] == "exon":
									if not locus in exon:
										exon[locus] = []
									exon[locus].append((fields[3],fields[4],fields[6]))
									if "frameshifts" in fields[8]:
										warnings.append("Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])

						newList = sorted(exon[locus], key=lambda x: x[1])
						exon[locus] = newList


						#Reconstruct CDS  ***********************************************************
						cdsSeq = ""
						if exon[locus][0][2]=="+":   #************* Positive strand
							for item in exon[locus]:
								cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
							cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

						else: # *********************** Negative Strand
							for a in range(len(exon[locus])-1,-1,-1):
								cdsSeq+= reverseComplement(genomeSeq[ int(exon[locus][a][0])-1: int(exon[locus][a][1]) ])
							cdsSeq += reverseComplement(genomeSeq[int(exon[locus][a][0])-4 :int(exon[locus][a][0])-1])


						notes = ""

						# Check CDS integrity *********************************************************
						foundStartCodon = True
						foundStopCodon = True
						if not  cdsSeq[:3]=="ATG" or not (cdsSeq[-3:]=="TGA" or cdsSeq[-3:]=="TAA" or cdsSeq[-3:]=="TAG" ):

							notes = "\n"+locus+"\n"
							notes += "Either the start or the stop codon was not found. Searching nearby....\n"




							# Look for ATG at the beginning of the sequence or closely ********************
							if exon[locus][0][2]=="+":   #************* Positive strand
								if not cdsSeq[:3]=="ATG":
									foundStartCodon = False
									print("- Looking for start codon upstream....\n")
									for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
										newStart = genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
										if newStart == "ATG":
											notes += "Valid start codon found upstream!\n"
											exon[locus][0]=(int(exon[locus][0][0])-a*3-3,exon[locus][0][1],exon[locus][0][2])
											gene[locus] = (int(exon[locus][0][0])-a*3-3,int(gene[locus][1]),gene[locus][2])
											cdsSeq = ""
											for item in exon[locus]:
												cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
											foundStartCodon = True
											notes += "- Start codon refined  "+str(a)+" codons upstream\n"
											numCodonRefines = a
											break
										if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
											foundStartCodon = False
											notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
											break
									#If the new start codon was not found in the region upstream then the downstream region is searched
									if foundStartCodon == False:
										notes += "- Looking for start codon downstream....\n"
										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											newStart = genomeSeq[int(exon[locus][0][0])-1+a*3:int(exon[locus][0][0])-1+a*3+3]
											if newStart == "ATG":
												notes += "Valid start codon found downstream!\n"
												exon[locus][0]=(int(exon[locus][0][0])+a*3,exon[locus][0][1],exon[locus][0][2])
												gene[locus] = (int(exon[locus][0][0])+a*3,int(gene[locus][1]),gene[locus][2])
												cdsSeq = ""
												for item in exon[locus]:
													cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
												foundStartCodon = True
												notes += "- Start codon refined  "+str(a)+" codons downstream\n"
												numCodonRefines = a
												break
											if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
												foundStartCodon = False
												notes += "- Found a stop codon while searching for start codon downstream!\nStart codon could not be found downstream\n"
												break
										if foundStartCodon == False:
											notes += "Start codon could not be found at this stage\n"

							else: # *********************** Negative Strand
								if not cdsSeq[:3]=="ATG" or not (cdsSeq[:3] =="TTG" and locus=="RL6"): #RL6 start with alternative start codon
									foundStartCodon = False
									notes += "- Looking for start codon upstream....\n"
									for a in range(len(sequence)-int(int(len(cdsSeq)/3))+30):
										#print "New start codons"
										newStart = reverseComplement(genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
										#print newStart
										if newStart == "ATG":
											notes += "Valid start codon found upstream!\n"
											exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])+a*3+3,exon[locus][-1][2])
											gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])+a*3+3,gene[locus][2])
											cdsSeq = ""
											for a1 in range(len(exon[locus])-1,-1,-1):
												cdsSeq+=reverseComplement(genomeSeq[ int(exon[locus][a1][0])-1: int(exon[locus][a1][1]) ])
											foundStartCodon = True
											notes += "- Start codon refined  "+str(a)+" codons upstream\n"
											numCodonRefines = a
											break
										if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
											foundStartCodon = False
											notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
											break
								#If the new start codon was not found in the region upstream then the downstream region is searched
									if foundStartCodon == False:
										notes += "- Looking for start codon downstream....\n"

										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											#print "New start codons"
											newStart = reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3-3:int(exon[locus][-1][1])-a*3])

											#print newStart
											if newStart == "ATG":
												notes += "Valid start codon found downstream!\n"
												exon[locus][-1]=(int(exon[locus][-1][0]),int(exon[locus][-1][1])-a*3,exon[locus][-1][2])
												gene[locus] = (int(gene[locus][0]),int(exon[locus][-1][1])-a*3,gene[locus][2])
												cdsSeq = ""
												for a1 in range(len(exon[locus])-1,-1,-1):
													cdsSeq+=reverseComplement(genomeSeq[ int(exon[locus][a1][0])-1: int(exon[locus][a1][1]) ])
												foundStartCodon = True
												notes += "- Start codon refined  "+str(a)+" codons downstream\n"
												numCodonRefines = a
												break
										if newStart == "TGA" or newStart == "TAA" or newStart=="TAG":
											foundStartCodon = False
											notes += "- Found a stop codon while searching for start codon downstream! \nStart codon could not be found downstream\n"
											break

										if foundStartCodon == False:
											notes += "Start codon could not be found at this stage\n"




							# Look for Stop codon at the end of the sequence or closely ********************
							if exon[locus][0][2]=="+":   #************* Positive strand
								#print "We start from", cdsSeq[-3:]
								if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
									notes += "- Looking for stop codon downstream....\n"
									foundStopCodon = False
									for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
										newStop = genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
										if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
											notes += "Valid stop codon found downstream\n"
											exon[locus][0]=(int(exon[locus][0][0]),int(exon[locus][-1][1])+a*3,exon[locus][0][2])
											gene[locus] = (int(gene[locus][0]), int(exon[locus][-1][1])+a*3, gene[locus][2])
											cdsSeq = ""
											for item in exon[locus]:
												cdsSeq+=genomeSeq[int(item[0])-1:int(item[1])]
											foundStopCodon = True
											notes += "- Stop codon refined " + \
												str(a)+" codon downstream\n"
											break


							else: # *********************** Negative Strand
								if not cdsSeq[-3:]=="TGA" and not cdsSeq[-3:]=="TAA" and not cdsSeq[-3:]=="TAG":
									notes += "- Looking for stop codon downstream....\n"
									for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
										#print "New Stop codons"
										newStop = reverseComplement(genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])
										#print a,newStop
										if newStop == "TAA" or newStop=="TGA" or newStop=="TAG":
											notes += "Valid stop codon found downstream\n"
											exon[locus][0]=(int(exon[locus][0][0])-a*3,exon[locus][0][1],exon[locus][0][2])
											gene[locus] =(int(exon[locus][0][0])-a*3, int(gene[locus][1]),gene[locus][2])
											cdsSeq = ""
											for b in range(len(exon[locus])-1,-1,-1):
												cdsSeq+=reverseComplement(genomeSeq[ int(exon[locus][b][0])-1: int(exon[locus][b][1]) ])
											foundStopCodon = True
											notes += "- Stop codon refined " + \
												str(a)+" codon(s) downstream\n"
											break



						if foundStartCodon == True and foundStopCodon == True and numCodonRefines <5 and abs(int(len(cdsSeq)/3) - len(sequence)) <10: # Write gff and cds file ********************
							#warnFile.write("A valid ORF for gene "+locus+" after prediction refinement\n")
							if  notes == "":
								notes = "A valid ORF has been found for gene "+locus+"!\n"
							else:
								notes += "A valid ORF has been found for gene "+locus+"!\n"
							gffNote = ""
							#  ******************* Check CDS integrity
							cdsGood = True
							plausablePrediction = True
							for a in range(0,len(cdsSeq)-3,+3):
								if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
									# Check if the shorter sequence is compatible with one of the models
									newProtLen = float(a)/3.0
									plausablePrediction = False
									newProt = Seq(cdsSeq[:a+3]).translate()
									for protein in protSeqs:
										score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
										if score / float(len(protSeqs[protein].seq)) >= 0.8:
											plausablePrediction = True
											if exon[locus][0][2]=="+":  #Check it if the strand is positive
												newmRNALength = 0
												newExonSet = {}
												if not locus in newExonSet:
													newExonSet[locus] = []
												for item in exon[locus]:
													if int(item[1])-int(item[0]) + newmRNALength > a+3:
														newExonSet[locus].append((int(item[0]),int(item[0]) + a - newmRNALength -1 ,item[2]))
														exon[locus] = newExonSet[locus]
														gene[locus] = (int(gene[locus][0]),int(item[0]) + a - newmRNALength -1 ,gene[locus][2])
														break
													else:
														newmRNALength += int(item[1]) - int(item[0])
														newExonSet[locus].append((int(item[0]),int(item[1]),item[2]))
												cdsSeq = cdsSeq[:a+3]
											else:  #Check it if the strand is negative
												newmRNALength = 0
												newExonSet = {}
												if not locus in newExonSet:
													newExonSet[locus] = []
												for a in range(len(exon[locus])-1,-1,-1):
													if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
														newExonSet[locus].append((int(exon[locus][a][1]) - a  + newmRNALength, int(exon[locus][a][1]),exon[locus][a][2]))
														#print "Previous exon locus",exon[locus]
														exon[locus] = newExonSet[locus]
														#print "after exon locus",exon[locus]
														#print gene[locus]
														#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
														break

													else:
														newmRNALength = int(exon[locus][a][1]) - int(exon[locus][a][0])
														newExonSet[locus].append((int(exon[locus][a][0]),int(exon[locus][a][1]),exon[locus][a][2]))
												cdsSeq = cdsSeq[:a+3]
												#print exon[locus][0][1]
												exon[locus][0] = (exon[locus][0][0], exon[locus][0][1]-6,exon[locus][0][2])
												#print exon[locus][0][1]
											break
										else:
											plausablePrediction = False

									if plausablePrediction == False:
										break


							if plausablePrediction == False:
								cdsGood = False
								gffNote = "note=Stop codon interrupts coding sequence. "
								notes += "- The coding sequence is interrupted by a stop codon\n"


							if not len(cdsSeq)%3 == 0:
								cdsGood = False
								notes += "- The coding sequence length is not multiple of 3\n"
								if not gffNote == "" :
									gffNote = "note= The coding sequence length is not multiple of 3"
								else:
									gffNote +=". The coding sequence length is not multiple of 3"



							if cdsGood == True:  # CDS passed quality check
								if exon[locus][0][2]=="+":   #************* Positive strand
									gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
									gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
									numExon = 1
									for item in exon[locus]:
										if item == exon[locus][-1]: #If this is the last exon include the stop codon in the coordinates
											gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(int(item[1])+3)+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
										else:
											gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
										numExon += 1


								else: # *********************** Negative Strand
									gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][-1][1]))+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
									gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][-1][1]))+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
									numExon = 1
									for item in exon[locus]:
										if item == exon[locus][0]:# and len(exon[locus]) == 1: #If this is the first exon include the stop codon in the coordinates
											gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(int(item[0])-3)+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
										else:
											gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
										numExon += 1

								cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

							else: # CDS DID NOT passed quality check
								#warnFile.write(notes+"\n")
								gffNote += "Pseudo "
								if exon[locus][0][2]=="+":   #************* Positive strand
									gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\t"+gffNote+"-gene\n")
									gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
									numExon = 1
									for item in exon[locus]:
										gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
										numExon += 1


								else: # *********************** Negative Strand
									gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
									gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
									numExon = 1
									for item in exon[locus]:
										gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
										numExon += 1

								cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

							warnFile.write(notes+"\n")




						else:
						#if foundStartCodon == False or foundStopCodon == False or abs(int(len(cdsSeq)/3) - len(sequence))>10:
							notes += "The prediction was not successfull! Now attempting a refinement....\n"
							self.logArea.configure(state='normal')
							self.logArea.insert(tk.END, "Annotation needs refinement....\n")
							self.logArea.see(tk.END)
							self.logArea.configure(state='disabled')
							self.logArea.update()
							print("Annotation needs refinement")
							exResult.close()
							#  ***************************************************************************
							#  ************************* Annotation refinement ***************************
							#  ***************************************************************************


							os.system(installationDirectory+"src/conda/bin/exonerate --model protein2genome tempFasta.fasta "+genomeName+" --showtargetgff -s 0 -n 1 --refine full --forcegtag --minintron 35 --maxintron 10000 >outputExonerate")

							#Check Exonerate output *****************************************************************
							#Check the proteins gave a match in the target genome
							exResult = open("outputExonerate")
							line = exResult.readline().rstrip()
							while not "Query range:" in line:
								line = exResult.readline().rstrip()
								if line is None:
									print("WARNING no protein found for ", locus)
									warnings.append(
										"Missing gene: locus "+locus+" did not provide any alignment")
							exResult.close()

							#Reconstruct Exons  **********************************************************************
							exResult = open("outputExonerate")
							while not line == "# --- START OF GFF DUMP ---":
								line = exResult.readline().rstrip()
								if line is None:
									#print "WARNING no protein found for ", locus
									warnings.append(
										"Missing gene: locus "+locus+" did not provide any alignment")
							for a in range(10):
								line = exResult.readline().rstrip()

							gene = {}
							exon = {}
							warnings = []
							#Collect the exonerate output
							while not line == "# --- END OF GFF DUMP ---":
								line = exResult.readline().rstrip()
								if line is None:
									#print "WARNING no protein found for ", locus
									warnings.append(
										"Missing gene: locus "+locus+" did not provide any alignment")
								fields = line.split("\t")
								if not line == "# --- END OF GFF DUMP ---":
									if fields[2] == "gene":
										if not locus in gene:
											gene[locus] = (
												fields[3], fields[4], fields[6])
										else:
											if (int(fields[4]) - int(fields[3])) > (int(int(gene[locus][1])) - int(int(gene[locus][0]))):
												gene[locus] = (
													fields[3], fields[4], fields[6])

									if fields[2] == "exon":
										if not locus in exon:
											exon[locus] = []
										exon[locus].append(
											(fields[3], fields[4], fields[6]))
										if "frameshifts" in fields[8]:
											warnings.append(
												"Frameshifts in exon "+str(fields[3])+" "+str(fields[4])+" "+fields[8])

							newList = sorted(exon[locus], key=lambda x: x[1])
							exon[locus] = newList

							#Reconstruct CDS  ***********************************************************
							cdsSeq = ""
							if exon[locus][0][2] == "+":  # ************* Positive strand
								for item in exon[locus]:
									cdsSeq += genomeSeq[int(item[0]) -
														1:int(item[1])]
								cdsSeq += genomeSeq[int(item[1]):int(item[1])+3]

							else:  # *********************** Negative Strand
								for a in range(len(exon[locus])-1, -1, -1):
									cdsSeq += reverseComplement(
										genomeSeq[int(exon[locus][a][0])-1: int(exon[locus][a][1])])
								cdsSeq += reverseComplement(
									genomeSeq[int(exon[locus][a][0])-4:int(exon[locus][a][0])-1])



							# Check CDS integrity *********************************************************
							foundStartCodon = True
							foundStopCodon = True
							if not cdsSeq[:3] == "ATG" or not (cdsSeq[-3:] == "TGA" or cdsSeq[-3:] == "TAA" or cdsSeq[-3:] == "TAG"):


								notes += "Either the start or the stop codon was not found. Searching nearby....\n"

								# Look for ATG at the beginning of the sequence or closely ********************
								if exon[locus][0][2] == "+":  # ************* Positive strand
									if not cdsSeq[:3] == "ATG":
										foundStartCodon = False
										print("- Looking for start codon upstream....\n")
										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											newStart = genomeSeq[int(
												exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3]
											if newStart == "ATG":
												notes += "Valid start codon found upstream!\n"
												exon[locus][0] = (
													int(exon[locus][0][0])-a*3-3, exon[locus][0][1], exon[locus][0][2])
												gene[locus] = (
													int(exon[locus][0][0])-a*3-3, int(gene[locus][1]), gene[locus][2])
												cdsSeq = ""
												for item in exon[locus]:
													cdsSeq += genomeSeq[int(
														item[0])-1:int(item[1])]
												foundStartCodon = True
												notes += "- Start codon refined  " + \
													str(a)+" codons upstream\n"
												numCodonRefines = a
												break
											if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
												foundStartCodon = False
												notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
												break
										#If the new start codon was not found in the region upstream then the downstream region is searched
										if foundStartCodon == False:
											notes += "- Looking for start codon downstream....\n"
											for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
												newStart = genomeSeq[int(
													exon[locus][0][0])-1+a*3+3:int(exon[locus][0][0])-1+a*3]
												if newStart == "ATG":
													notes += "Valid start codon found downstream!\n"
													exon[locus][0] = (
														int(exon[locus][0][0])+a*3+3, exon[locus][0][1], exon[locus][0][2])
													gene[locus] = (
														int(exon[locus][0][0])+a*3+3, int(gene[locus][1]), gene[locus][2])
													cdsSeq = ""
													for item in exon[locus]:
														cdsSeq += genomeSeq[int(
															item[0])-1:int(item[1])]
													foundStartCodon = True
													notes += "- Start codon refined  " + \
														str(a) + \
														" codons upstream\n"
													numCodonRefines = a
													break
												if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
													foundStartCodon = False
													notes += "- Found a stop codon while searching for start codon downstream!\nStart codon could not be found downstream\n"
													break

								else:  # *********************** Negative Strand
									# RL6 start with alternative start codon
									if not cdsSeq[:3] == "ATG" or not (cdsSeq[:3] == "TTG" and locus == "RL6"):
										foundStartCodon = False
										notes += "- Looking for start codon upstream....\n"
										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											#print "New start codons"
											newStart = reverseComplement(
												genomeSeq[int(exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3])
											#print newStart
											if newStart == "ATG":
												notes += "Valid start codon found upstrean!\n"
												exon[locus][-1] = (int(exon[locus][-1][0]), int(
													exon[locus][-1][1])+a*3+3, exon[locus][-1][2])
												gene[locus] = (int(gene[locus][0]), int(
													exon[locus][-1][1])+a*3+3, gene[locus][2])
												cdsSeq = ""
												for a1 in range(len(exon[locus])-1, -1, -1):
													cdsSeq += reverseComplement(
														genomeSeq[int(exon[locus][a1][0])-1: int(exon[locus][a1][1])])
												foundStartCodon = True
												notes += "- Start codon refined  " + \
													str(a)+" codons upstream\n"
												numCodonRefines = a
												break
											if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
												foundStartCodon = False
												notes += "- Found a stop codon while searching for start codon upstream!\nStart codon could not be found upstream\n"
												break
									#If the new start codon was not found in the region upstream then the downstream region is searched
										if foundStartCodon == False:
											notes += "- Looking for start codon downstream....\n"
											for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
												#print "New start codons"
												newStart = reverseComplement(genomeSeq[int(exon[locus][-1][1])-a*3-3:int(exon[locus][-1][1])-a*3])
												#print newStart
												if newStart == "ATG":
													notes += "Valid start codon found dowstream!\n"
													exon[locus][-1] = (int(exon[locus][-1][0]), int(exon[locus][-1][1])-a*3, exon[locus][-1][2])
													gene[locus] = (int(gene[locus][0]), int(
														exon[locus][-1][1])-a*3-3, gene[locus][2])
													cdsSeq = ""
													for a1 in range(len(exon[locus])-1, -1, -1):
														cdsSeq += reverseComplement(
															genomeSeq[int(exon[locus][a1][0])-1: int(exon[locus][a1][1])])
													foundStartCodon = True
													notes += "- Start codon refined  " + \
														str(a) + \
														" codons dowstream\n"
													numCodonRefines = a
													break
												if newStart == "TGA" or newStart == "TAA" or newStart == "TAG":
													foundStartCodon = False
													notes += "- Found a stop codon while searching for start codon downstream! \nStart codon could not be found downstream\n"
													break
										if foundStartCodon == False:
											notes += "Start codon could not be found at this stage\n"


								# Look for Stop codon at the end of the sequence or closely ********************
								if exon[locus][0][2] == "+":  # ************* Positive strand
									if not cdsSeq[-3:] == "TGA" and not cdsSeq[-3:] == "TAA" and not cdsSeq[-3:] == "TAG":
										notes += "- Looking for stop codon downstream....\n"
										foundStopCodon = False
										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											newStop = genomeSeq[int(
												exon[locus][-1][1])+a*3:int(exon[locus][-1][1])+a*3+3]
											if newStop == "TAA" or newStop == "TGA" or newStop == "TAG":
												notes += "Valid stop codon found downstream!\n"
												exon[locus][0] = (int(exon[locus][0][0]), int(
													exon[locus][-1][1])+a*3, exon[locus][0][2])
												gene[locus] = (int(gene[locus][0]), int(
													exon[locus][-1][1])+a*3, gene[locus][2])
												cdsSeq = ""
												for item in exon[locus]:
													cdsSeq += genomeSeq[int(
														item[0])-1:int(item[1])]
												foundStopCodon = True
												notes += "- Stop codon refined " + \
													str(a) + \
													" codon downstream\n"
												break

								else:  # *********************** Negative Strand
									if not cdsSeq[-3:] == "TGA" and not cdsSeq[-3:] == "TAA" and not cdsSeq[-3:] == "TAG":
										foundStopCodon = False
										notes += "- Looking for stop codon downstream....\n"
										for a in range(len(sequence)-int(len(cdsSeq)/3)+30):
											#print "New Stop codons"
											newStop = reverseComplement(
												genomeSeq[int(exon[locus][0][0])-1-a*3-3:int(exon[locus][0][0])-1-a*3])

											if newStop == "TAA" or newStop == "TGA" or newStop == "TAG":
												notes += "Valid stop codon found downstream!\n"

												exon[locus][0] = (
													int(exon[locus][0][0])-a*3, exon[locus][0][1], exon[locus][0][2])
												gene[locus] = (
													int(exon[locus][0][0])-a*3, int(gene[locus][1]), gene[locus][2])

												cdsSeq = ""
												for b in range(len(exon[locus])-1, -1, -1):
													cdsSeq += reverseComplement(
														genomeSeq[int(exon[locus][b][0])-1: int(exon[locus][b][1])])
												foundStopCodon = True
												notes += "- Stop codon refined " + \
													str(a) + \
													" codon downstream\n"
												break

								warnFile.write(notes)#+"\n")

							# Write gff and cds file ********************
							if foundStartCodon == True and foundStopCodon == True:

								warnFile.write("A valid ORF for gene "+locus+" after prediction refinement!\n")
								if notes == "":
									notes = "A valid ORF for gene "+locus+" after prediction refinement!\n"
								else:
									notes += "A valid ORF has been found for gene "+locus+" after prediction refinement!\n"
								gffNote = ""
								#  ******************* Check CDS integrity
								cdsGood = True
								plausablePrediction = True

								for a in range(0, len(cdsSeq)-3, +3):
									plausablePrediction = True
									if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
										# Check if the shorter sequence is compatible with one of the models
										newProtLen = float(a)/3.0
										plausablePrediction = False
										newProt = Seq(cdsSeq[:a+3]).translate()
										for protein in protSeqs:
											score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
											if score / float(len(protSeqs[protein].seq)) >= 0.8:
												plausablePrediction = True
												# Check it if the strand is positive
												if exon[locus][0][2] == "+":
													newmRNALength = 0
													newExonSet = {}
													if not locus in newExonSet:
														newExonSet[locus] = []
													for item in exon[locus]:
														if int(item[1])-int(item[0]) + newmRNALength > a+3:
															newExonSet[locus].append(
																(int(item[0]), int(item[0]) + a - newmRNALength - 1, item[2]))
															exon[locus] = newExonSet[locus]
															gene[locus] = (int(gene[locus][0]), int(
																item[0]) + a - newmRNALength - 1, gene[locus][2])
															break
														else:
															newmRNALength += int(
																item[1]) - int(item[0])
															newExonSet[locus].append(
																(int(item[0]), int(item[1]), item[2]))
													cdsSeq = cdsSeq[:a+3]
												else:  # Check it if the strand is negative
													newmRNALength = 0
													newExonSet = {}
													if not locus in newExonSet:
														newExonSet[locus] = []
													for a in range(len(exon[locus])-1, -1, -1):
														if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
															newExonSet[locus].append(
																(int(exon[locus][a][1]) - a + newmRNALength, int(exon[locus][a][1]), exon[locus][a][2]))
															#print "Previous exon locus",exon[locus]
															exon[locus] = newExonSet[locus]
															#print "after exon locus",exon[locus]
															#print gene[locus]
															#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
															break

														else:
															newmRNALength = int(
																exon[locus][a][1]) - int(exon[locus][a][0])
															newExonSet[locus].append(
																(int(exon[locus][a][0]), int(exon[locus][a][1]), exon[locus][a][2]))
													cdsSeq = cdsSeq[:a+3]
													#print exon[locus][0][1]
													exon[locus][0] = (
														exon[locus][0][0], exon[locus][0][1]-6, exon[locus][0][2])
													#print exon[locus][0][1]
												break
											else:
												plausablePrediction = False

										if plausablePrediction == False:
											break




								if plausablePrediction==False:

									cdsGood = False
									gffNote = "note=Stop codon interrupts coding sequence. "
									notes += "- The coding sequence is interrupted by a stop codon\n"
									#break

								if not len(cdsSeq) % 3 == 0:
									cdsGood = False
									notes += "- The coding sequence length is not multiple of 3\n"
									if not gffNote == "":
										gffNote = "note= The coding sequence length is not multiple of 3"
									else:
										gffNote += ". The coding sequence length is not multiple of 3"

								if cdsGood == True:  # CDS passed quality check
									if exon[locus][0][2]=="+":   #************* Positive strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
										numExon = 1
										for item in exon[locus]:
											if item == exon[locus][-1]: #If this is the last exon include the stop codon in the coordinates
												gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(int(item[1])+3)+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											else:
												gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											numExon += 1


									else: # *********************** Negative Strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][-1][1])+3)+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][-1][1])+3)+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
										numExon = 1
										for item in exon[locus]:
											if item == exon[locus][0]:# and len(exon[locus]) == 1: #If this is the first exon include the stop codon in the coordinates
												gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(int(item[0])-3)+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											else:
												gffFile.write(assemblyName+"\texonerate\t"+"CDS"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											numExon += 1

									cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

								else: # CDS DID NOT passed quality check

									warnFile.write(notes+"\n")
									gffNote += "Pseudo "
									if exon[locus][0][2]=="+":   #************* Positive strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\t"+gffNote+"-gene\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
										numExon = 1
										for item in exon[locus]:
											gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
											numExon += 1


									else: # *********************** Negative Strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
										numExon = 1
										for item in exon[locus]:
											gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
											numExon += 1



									cdsFile.write(">"+locus+" " +gffNote+"-CDS\n"+cdsSeq+"\n")

							else:
								warnFile.write("An incomplete ORF for gene "+locus+" after prediction refinement!\n")
								if notes == "":
									notes = "An incomplete ORF for gene "+locus+" after prediction refinement!\n"
								else:
									notes += "An imcoplete ORF has been found for gene "+locus+" after prediction refinement!\n"
								gffNote = ""
								#  ******************* Check CDS integrity
								cdsGood = True
								plausablePrediction = True
								for a in range(0, len(cdsSeq)-3, +3):
									if cdsSeq[a:a+3] == "TAA" or cdsSeq[a:a+3] == "TGA" or cdsSeq[a:a+3] == "TAG":
										# Check if the shorter sequence is compatible with one of the models
										newProtLen = float(a)/3.0
										plausablePrediction = False
										newProt = Seq(cdsSeq[:a+3]).translate()
										for protein in protSeqs:
											score = pairwise2.align.localxx(newProt, protSeqs[protein].seq,score_only=True)
											if score / float(len(protSeqs[protein].seq)) >= 0.8:
												plausablePrediction = True
												# Check it if the strand is positive
												if exon[locus][0][2] == "+":
													newmRNALength = 0
													newExonSet = {}
													if not locus in newExonSet:
														newExonSet[locus] = []
													for item in exon[locus]:
														if int(item[1])-int(item[0]) + newmRNALength > a+3:
															newExonSet[locus].append(
																(int(item[0]), int(item[0]) + a - newmRNALength - 1, item[2]))
															exon[locus] = newExonSet[locus]
															gene[locus] = (int(gene[locus][0]), int(
																item[0]) + a - newmRNALength - 1, gene[locus][2])
															break
														else:
															newmRNALength += int(
																item[1]) - int(item[0])
															newExonSet[locus].append(
																(int(item[0]), int(item[1]), item[2]))
													cdsSeq = cdsSeq[:a+3]
												else:  # Check it if the strand is negative
													newmRNALength = 0
													newExonSet = {}
													if not locus in newExonSet:
														newExonSet[locus] = []
													for a in range(len(exon[locus])-1, -1, -1):
														if int(exon[locus][a][1])-int(exon[locus][a][0]) + newmRNALength > a+3:
															newExonSet[locus].append(
																(int(exon[locus][a][1]) - a + newmRNALength, int(exon[locus][a][1]), exon[locus][a][2]))
															#print "Previous exon locus",exon[locus]
															exon[locus] = newExonSet[locus]
															#print "after exon locus",exon[locus]
															#print gene[locus]
															#gene[locus] = (int(exon[locus][a][1]) - a  + newmRNALength, int(gene[locus][1]), gene[locus][2])
															break

														else:
															newmRNALength = int(
																exon[locus][a][1]) - int(exon[locus][a][0])
															newExonSet[locus].append(
																(int(exon[locus][a][0]), int(exon[locus][a][1]), exon[locus][a][2]))
													cdsSeq = cdsSeq[:a+3]
													#print exon[locus][0][1]
													exon[locus][0] = (
														exon[locus][0][0], exon[locus][0][1]-6, exon[locus][0][2])
													#print exon[locus][0][1]
												break
								if plausablePrediction==False:
									cdsGood = False
									gffNote = "note=Stop codon interrupts coding sequence. "
									notes += "- The coding sequence is interrupted by a stop codon\n"


								if not len(cdsSeq) % 3 == 0:
									cdsGood = False
									notes += "- The coding sequence length is not multiple of 3\n"
									if not gffNote == "":
										gffNote = "note= The coding sequence length is not multiple of 3"
									else:
										gffNote += ". The coding sequence length is not multiple of 3"

								if cdsGood == True:  # CDS passed quality check
									if exon[locus][0][2]=="+":   #************* Positive strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][0][0]))+"\t"+str(int(exon[locus][-1][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
										numExon = 1
										for item in exon[locus]:
											if item == exon[locus][-1]: #If this is the last exon include the stop codon in the coordinates
												gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(int(item[1])+3)+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											else:
												gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											numExon += 1


									else: # *********************** Negative Strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(exon[locus][-1][1])+3)+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(exon[locus][-1][1])+3)+"\t"+str(int(exon[locus][0][0]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+"\n")
										numExon = 1
										for item in exon[locus]:
											if item == exon[locus][0]:# and len(exon[locus]) == 1: #If this is the first exon include the stop codon in the coordinates
												gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(int(item[0])-3)+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											else:
												gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+"\n")
											numExon += 1

									cdsFile.write(">"+locus+" +\n"+cdsSeq+"\n")

								else: # CDS DID NOT passed quality check
									warnFile.write(notes+"\n")
									gffNote += "Pseudo "
									if exon[locus][0][2]=="+":   #************* Positive strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+"\t"+gffNote+"-gene\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0]))+"\t"+str(int(gene[locus][1])+3)+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
										numExon = 1
										for item in exon[locus]:
											gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
											numExon += 1


									else: # *********************** Negative Strand
										gffFile.write(assemblyName+"\texonerate\t"+"gene"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_gene;Name="+locus+";Product="+locus+";"+gffNote+"-gene\n")
										gffFile.write(assemblyName+"\texonerate\t"+"mRNA"+"\t"+str(int(gene[locus][0])-3)+"\t"+str(int(gene[locus][1]))+"\t.\t"+str(gene[locus][2])+"\t.\tID="+locus+"_mRNA;Parent="+locus+"_gene;Name="+locus+".1;Product="+locus+";"+gffNote+"-mRNA\n")
										numExon = 1
										for item in exon[locus]:
											gffFile.write(assemblyName+"\texonerate\t"+"misc_feature"+"\t"+str(item[0])+"\t"+str(item[1])+"\t.\t"+str(item[2])+"\t.\tID="+locus+"_cds"+str(numExon)+";Parent="+locus+"_mRNA;Name="+locus+".1;Product="+locus+";"+gffNote+"-CDS\n")
											numExon += 1



									cdsFile.write(">"+locus+" " +gffNote+"-CDS\n"+cdsSeq+"\n")













					os.system("rm tempFasta.fasta outputExonerate outputBlast.txt -f")
					cdsFile.close()
					gffFile.close()
					warnFile.close()

					#Translate valid cds in proteins

					sequences = SeqIO.to_dict(SeqIO.parse(suffixName+"_cds.fasta","fasta"))

					for seq in sequences:
						if not "pseudo" in str(sequences[seq].description):
							protSeq = (sequences[seq].seq).translate()
							protFile.write(">"+str(sequences[seq].description)+"\n"+str(protSeq)+"\n")

					protFile.close()


					




					##refine gff file
					geneInfo = {}
					mRNAInfo = {}
					cdsInfo = {}
					cdsExons = {}
					geneCoord = {}
					geneDirection = {}
					exonCoordAll = {}

					inputFileGFF = suffixName+"_annotation.gff"

					infileGFF = open(inputFileGFF)
					line = infileGFF.readline()
					fields = line.split("\t")
					genomeName = fields[0]
					infileGFF.close()

					infileGFF = open(inputFileGFF)
					while True:
						line = infileGFF.readline().rstrip()
						if not line:
							break
						fields = line.split("\t")


						if fields[2] == "gene":
							infoField = fields[-1]
							if  infoField[:4]=="note":
								infoField = fields[-2]
							locus = ((infoField.split("ID="))[-1].split("_"))[0]
							if not locus in geneInfo:
								geneInfo[locus] = infoField

							if not locus in geneDirection:
								geneDirection[locus]  = (fields[5],fields[6],fields[7])



						if fields[2] == "mRNA":
							infoField = fields[-1]
							if infoField[:4]=="note":
								infoField = fields[-2]
							locus = ((infoField.split("ID="))[-1].split("_"))[0]
							if not locus in mRNAInfo:
								mRNAInfo[locus] = infoField


						if fields[2] == "CDS":
							infoField = fields[-1]
							if infoField[:4]=="note":
								infoField = fields[-2]
							locus = ((infoField.split("ID="))[-1].split("_"))[0]

							if not locus in cdsInfo:
								cdsInfo[locus] = []
							if int(fields[3]) > int(fields[4]):
								temp =  fields[3]
								fields[3] = fields[4]
								fields[4] = temp

							if not locus in exonCoordAll:
								exonCoordAll[locus] = []
							exonCoordAll[locus].append(int(fields[3]))
							exonCoordAll[locus].append(int(fields[4]))

							newCDSLine = ""
							for a in fields:
								newCDSLine +=a+"\t"

							cdsInfo[locus].append(newCDSLine)



						if fields[2] == "misc_feature":
							infoField = fields[-1]
							if infoField[:4]=="note":
								infoField = fields[-2]
							locus = ((infoField.split("ID="))[-1].split("_"))[0]


							if not locus in cdsInfo:
								cdsInfo[locus] = []
							if int(fields[3]) > int(fields[4]):
								temp =  fields[3]
								fields[3] = fields[4]
								fields[4] = temp

							if not locus in exonCoordAll:
								exonCoordAll[locus] = []
							exonCoordAll[locus].append(int(fields[3]))
							exonCoordAll[locus].append(int(fields[4]))

							newCDSLine = ""
							for a in fields:
								newCDSLine +=a+"\t"

							cdsInfo[locus].append(newCDSLine)





					for cds in exonCoordAll:
						geneCoord[cds] = (min(exonCoordAll[cds]),max(exonCoordAll[cds]))

					#rewrite the gfffile
					outfile = open("temp140875.gff","w")
					for locus in geneCoord:
						outfile.write(genomeName+"\tGRACy\tgene\t"+str(geneCoord[locus][0])+"\t"+str(geneCoord[locus][1])+"\t")
						for item in geneDirection[locus]:
							outfile.write(item+"\t")
						outfile.write(geneInfo[locus]+"\n")

						if not "misc_feature" in cdsInfo[locus][0]:
							outfile.write(genomeName+"\tGRACy\tmRNA\t"+str(geneCoord[locus][0])+"\t"+str(geneCoord[locus][1])+"\t")
							for item in geneDirection[locus]:
								outfile.write(item+"\t")
							outfile.write(mRNAInfo[locus]+"\n")

						for cds in cdsInfo[locus]:
							outfile.write(cds+"\n")

					

					outfile.close()
					os.system("mv temp140875.gff "+suffixName+"_annotation.gff")
					os.system("mv "+suffixName+"* "+outputFolder)
					infileGFF.close()





			self.logArea.configure(state='normal')
			self.logArea.insert(tk.END, "\n\nAnnotations completed!\n")
			self.logArea.see(tk.END)
			self.logArea.configure(state='disabled')
			self.logArea.update()




		'''This class configures and populates the toplevel window.
		   top is the toplevel containing window.'''
		_bgcolor = '#d9d9d9'  # X11 color: 'gray85'
		_fgcolor = '#000000'  # X11 color: 'black'
		_compcolor = '#d9d9d9' # X11 color: 'gray85'
		_ana1color = '#d9d9d9' # X11 color: 'gray85'
		_ana2color = '#ececec' # Closest X11 color: 'gray92'
		self.style = ttk.Style()
		if sys.platform == "win32":
			self.style.theme_use('winnative')
		self.style.configure('.',background=_bgcolor)
		self.style.configure('.',foreground=_fgcolor)
		self.style.configure('.',font="TkDefaultFont")
		self.style.map('.',background=
			[('selected', _compcolor), ('active',_ana2color)])

		w=900
		h=450
		ws = root.winfo_screenwidth() # width of the screen
		hs = root.winfo_screenheight() # height of the screen
		# calculate x and y coordinates for the Tk root window
		x = (ws/2) - (w/2)
		y = (hs/2) - (h/2)

		top.geometry('%dx%d+%d+%d' % (w, h, x, y))
		top.title("Annotation tool")
		top.configure(highlightcolor="black",background="white")

		self.inputFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.inputFileLabel.configure(text="Input file")
		self.inputFileLabel.place(x=20,y=20,width=70,height=20)

		self.inputFileEntry = tk.Entry(top)
		self.inputFileEntry.place(x=20,y=40,width=750,height=30)
		self.inputFileEntry.insert(0,"Please select a genomes list....")

		self.inputFileButton = tk.Button(top,command=openInputFile,background="#204949",foreground="white")
		self.inputFileButton.place(x=780,y=40,width=100,height=30)
		self.inputFileButton.configure(text="Open file")


		self.outputFolderLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.outputFolderLabel.configure(text="Output folder")
		self.outputFolderLabel.place(x=20,y=100,width=90,height=20)

		self.outputFolderEntry = tk.Entry(top)
		self.outputFolderEntry.place(x=20,y=120,width=750,height=30)
		self.outputFolderEntry.insert(0,"Please select an output folder....")

		self.outputFolderButton = tk.Button(top,command=openOutputFolder,background="#204949",foreground="white")
		self.outputFolderButton.place(x=780,y=120,width=100,height=30)
		self.outputFolderButton.configure(text="Open folder")



		self.runButton = tk.Button(top,command=mainAnnotationAlgorithm,background="#204949",foreground="white")
		self.runButton.place(x=720,y=400,width=160,height=30)
		self.runButton.configure(text="Run")

		#self.runButton = tk.Button(top,command=exitProgram,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		#self.runButton.place(x=780,y=350,width=100,height=30)
		#self.runButton.configure(text="Exit")

		self.logFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.logFileLabel.place(x=20,y=180,height=20,width=100)
		self.logFileLabel.configure(text="Log window")

		self.logFrame = tk.Frame(top)
		self.logFrame.place(x=20, y=200, height=230, width=670)
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(borderwidth="2")
		self.logArea = tk.Text(top,state='disabled')
		self.logArea.place(x=25,y=205,height=220, width=660)
		self.logArea.configure(background="white",borderwidth=5)
		self.logArea.configure(selectbackground="#c4c4c4")

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Annotation.jpg")
		photo1 =ImageTk.PhotoImage(image)
		logoLabel = tk.Label(top, compound=tk.TOP,height=200,width=200,image=photo1,borderwidth=0,highlightthickness=0)
		logoLabel.place(x=700,y=190)
		logoLabel.image = photo1




if __name__ == '__main__':
	vp_start_gui()
