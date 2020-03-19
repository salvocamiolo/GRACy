import sys
from PIL import ImageTk, Image
from matplotlib.patches import Patch
from matplotlib import colors as colorcode
from tkinter import simpledialog
from tkinter import font

installationDirectory = sys.argv[1]



from  PyPDF2 import PdfFileWriter, PdfFileReader
from io import StringIO
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import tkinter as tk
from tkinter import ttk


import tkinter
from tkinter import filedialog as tkFileDialog
import sys
import datetime
import random as rd
from Bio import SeqIO
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from Bio.SeqFeature import SeqFeature, FeatureLocation



def vp_start_gui():
	'''Starting point when module is the main routine.'''
	global val, w, root



	root = tk.Tk()
	top = Toplevel1 (root)




	#top.InputFileButton.bind('<Button-1>',openInputFile)

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
			self.InputFileEntry.delete(0,tk.END)
			self.InputFileEntry.insert(0,inputFile)

		def changeDBFile():
			dbfile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select a kmer database")
			dbfile = (dbfile.split("/"))[-1]
			self.dbEntry.delete(0,tk.END)
			self.dbEntry.insert(0,dbfile)

		def openOutputFolder():
			outputFolder = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
			self.OutputFolderEntry.delete(0,tk.END)
			self.OutputFolderEntry.insert(0,outputFolder)

		#**************************************************************************
		#***************** Main Genotyping   Algorithm Start **********************
		#**************************************************************************

		def runGenotyping():

			#Initialize variables
			orderedHyperLoci = ["rl5a","rl6","rl12","rl13","ul1","ul9","ul11","ul20","ul73","ul74","ul120","ul139","ul146"]

			#numReads = int(self.numReadsEntry.get())
			dbfile = installationDirectory+"src/scripts/genotyping/kmerDB/"+self.dbEntry.get()
			NumThreads = int(self.numThreadsEntry.get())
			inputFile = self.InputFileEntry.get()

			if not inputFile == "Please select file...." :
				if not self.OutputFolderEntry.get() =="Please select a folder....":
					os.system("mkdir "+self.OutputFolderEntry.get())
					datasets = []
					infile=open(inputFile)
					while True:
						read1 = infile.readline().rstrip()
						if not read1:
							break
						read2 = infile.readline().rstrip()
						if not read2:
							break

						datasets.append((read1,read2))

					numDatasets = len(datasets)
					progressBarIncrement = 100/(numDatasets*20)
					step = 0
					for dataset in datasets:
						read1 = dataset[0]
						read2 = dataset[1]

						#Extract path and files info
						inputPath =  '/'.join(read1.split('/')[0:-1])
						fileRoot1 = ((read1.split("/"))[-1])
						if fileRoot1[-2:] == "fq":
							fileRoot1 = fileRoot1[:-3]
						else:
							fileRoot1 = fileRoot1[:-6]

						fileRoot2 = ((read2.split("/"))[-1])
						if fileRoot2[-2:] == "fq":
							fileRoot2 = fileRoot2[:-3]
						else:
							fileRoot2 = fileRoot2[:-6]

						sampleRoot = fileRoot1[:-2]
						outfile = open(fileRoot1+"_IDCard.txt","w")
						logFile = open(fileRoot1+"_logFile.txt","w")
						now = datetime.datetime.now()
						logFile.write("Genotyping sample "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

						fileList = "list"+str(rd.randint(0,1000000))+".txt"
						listFile = open(fileList,"w")
						listFile.write(read1+"\n"+read2+"\n")
						listFile.close()

						#Perform deduplication
						now = datetime.datetime.now()
						logFile.write("Deduplication for dataset "+fileRoot1+"  started at "+now.strftime("%H:%M")+"\n")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "******************************************************************************\n")
						self.logArea.insert(tk.END, "Genotyping dataset "+fileRoot1+"/"+fileRoot2+"\n")
						self.logArea.insert(tk.END, "******************************************************************************\n\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Performing deduplication....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						#print "Performing deduplication...."
						dedupFile1 = fileRoot1+"_noDup_1.fq"
						dedupFile2 = fileRoot2+"_noDup_2.fq"
						os.system(installationDirectory+"src/conda/bin/fastuniq -i "+fileList+" -t q -o "+dedupFile1+" -p "+dedupFile2)
						os.system("rm -f "+fileList)

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						now = datetime.datetime.now()
						logFile.write("Deduplication for dataset "+sampleRoot+"  ended at "+now.strftime("%H:%M")+"\n\n")
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						#******************* Calculate average coverage for deduplicated reads   **************
						now = datetime.datetime.now()
						logFile.write("Reference coverage calculation for dataset "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Performing reference alignment....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/bowtie2 -1 "+dedupFile1+" -2 "+dedupFile2+" -x "+installationDirectory+"src/scripts/genotyping/fastaFiles/hcmvReference -S alignmenthsbfy43223.sam >null 2>&1")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Converting sam to bam....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/samtools view -bS -h alignmenthsbfy43223.sam > alignmenthsbfy43223.bam")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Sorting bam....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/samtools sort -o alignmenthsbfy43223_sorted.bam alignmenthsbfy43223.bam")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Calculating average coverage....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/samtools depth  alignmenthsbfy43223_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage.txt")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						avCovFile = open("avCoverage.txt")
						avCov = float(avCovFile.readline().rstrip())
						avCovFile.close()
						detectionTreshold = float(avCov*float(self.CutoffText.get()))
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Average coverage for deduplicated reads: "+str(avCov)+"\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system("rm -f alignmenthsbfy43223.* avCoverage.txt "+fileList)
						now = datetime.datetime.now()
						logFile.write("Reference coverage calculation for dataset "+sampleRoot+"  ended at "+now.strftime("%H:%M")+"\n\n")

						#Collect kmers from database
						geneKmers = {}
						kmerdbfile = open(dbfile)
						kmerdbfile.readline()
						while True:
							line = kmerdbfile.readline().rstrip()
							if not line:
								break
							fields = line.split("\t")
							if not (fields[0],fields[1]) in geneKmers:
								geneKmers[(fields[0],fields[1])] = []
								kmerseqs = fields[2].split(",")
								for item in kmerseqs:
									if not len(item) == 0:
										kmerLength = len(item)
										geneKmers[(fields[0],fields[1])].append(item)


						#Get sequences in memory
						reads = []
						numSeq = 0
						overallGeneInfo = {}
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Caclculate kmer frequencies for file "+dedupFile1+"....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile1+" -o "+dedupFile1+"_kmerCount.jf")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						numSeq = 0
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Caclculate kmer frequencies for file "+dedupFile2+"....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						os.system(installationDirectory+"src/conda/bin/jellyfish count -m 17 -s 100M -t 8 -C "+dedupFile2+" -o "+dedupFile2+"_kmerCount.jf")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Merging kmer files....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						os.system(installationDirectory+"src/conda/bin/jellyfish merge "+dedupFile1+"_kmerCount.jf "+dedupFile2+"_kmerCount.jf")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						step = step+1
						self.progressbar['value']=int( step*progressBarIncrement)
						self.progressbar.update()

						readsKmer = {}

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Loading kmers for sequences in "+dedupFile1+" into memory....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						for seq_record in SeqIO.parse(dedupFile1,"fastq"):
							sequence = str(seq_record.seq)
							revSequence = str(seq_record.seq.reverse_complement())
							reads.append(str(seq_record.seq))
							reads.append(str(seq_record.seq.reverse_complement()))
							for a in range(0,len(sequence)-kmerLength+1):
								if not sequence[a:a+kmerLength] in readsKmer:
									readsKmer[sequence[a:a+kmerLength]] = []
								readsKmer[sequence[a:a+kmerLength]].append(str(seq_record.id))

							for a in range(0,len(revSequence)-kmerLength+1):
								if not revSequence[a:a+kmerLength] in readsKmer:
									readsKmer[revSequence[a:a+kmerLength]] = []
								readsKmer[revSequence[a:a+kmerLength]].append(str(seq_record.id))

							numSeq +=1
							#if numSeq == int(numReads/2):
							#    break

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Loading kmers for sequences in "+dedupFile2+" into memory....")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						for seq_record in SeqIO.parse(dedupFile2,"fastq"):
							sequence = str(seq_record.seq)
							revSequence = str(seq_record.seq.reverse_complement())
							reads.append(str(seq_record.seq))
							reads.append(str(seq_record.seq.reverse_complement()))
							for a in range(0,len(sequence)-16):
								if not sequence[a:a+17] in readsKmer:
									readsKmer[sequence[a:a+17]] = []
								readsKmer[sequence[a:a+17]].append(str(seq_record.id))

							for a in range(0,len(revSequence)-16):
								if not revSequence[a:a+17] in readsKmer:
									readsKmer[revSequence[a:a+17]] = []
								readsKmer[revSequence[a:a+17]].append(str(seq_record.id))

							numSeq +=1
							#if numSeq == int(numReads/2):
							#    break

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Done!\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						now = datetime.datetime.now()
						logFile.write("Genotyping for dataset "+sampleRoot+"  started at "+now.strftime("%H:%M")+"\n")

						#Start searching for signatures
						for gene in orderedHyperLoci:
							self.logArea.configure(state='normal')
							self.logArea.insert(tk.END, "Genotyping gene "+gene+" for sample "+sampleRoot+"....")
							self.logArea.see(tk.END)
							self.logArea.configure(state='disabled')
							self.logArea.update()

							matchedReads = {}
							#Collect specific kmers for the genotypes of this gene
							specificKmerGroup = {}
							for item in geneKmers:
								if item[0] == gene:
									if not item[1] in specificKmerGroup:
										specificKmerGroup[item[1]] = []
									for seqs in geneKmers[item]:
										specificKmerGroup[item[1]].append(seqs)




							countSeq = {}
							totCount = 0
							numMatchedKmers = {}
							for gr in specificKmerGroup:


								if not gr in countSeq:
									countSeq[gr] = 0

								if not gr in matchedReads:
									matchedReads[gr] = set()

								if not gr in numMatchedKmers:
									numMatchedKmers[gr] = 0

								command = installationDirectory+"src/conda/bin/jellyfish query mer_counts_merged.jf "
								for querySeq in specificKmerGroup[gr]:
									command += querySeq
									command += " "
								command += " >counts.txt"

								os.system(command)

								countFile = open("counts.txt")


								while True:
									countLine = countFile.readline().rstrip()
									if not countLine:
										break
									countFields = countLine.split(" ")
									if int(countFields[1])>0:
										numMatchedKmers[gr] +=1
										for item in readsKmer[countFields[0]]:
											matchedReads[gr].add(item)
								countFile.close()
								logFile.write(gene+"\t"+"For genotype "+gr+" there are "+str(len(matchedReads[gr]))+" reads matching\n")
								if len(matchedReads[gr])>=detectionTreshold:
									countSeq[gr] = len(matchedReads[gr])
									totCount = totCount + len(matchedReads[gr])

							self.logArea.configure(state='normal')
							self.logArea.insert(tk.END, "Done!\n")
							self.logArea.see(tk.END)
							self.logArea.configure(state='disabled')
							self.logArea.update()

							outfile.write(gene)


							for gr in countSeq:
								if countSeq[gr]>0:
									percentage = float(countSeq[gr])/float(totCount)
									outfile.write("\t"+gr+"\t"+str(percentage)+"\t"+str(totCount)+"\t"+str(numMatchedKmers[gr]))
									self.logArea.configure(state='normal')
									self.logArea.insert(tk.END, "Found signature for genotype "+gr+" with a percentage of "+(str(percentage*100))[:5]+"%\n")
									self.logArea.see(tk.END)
									self.logArea.configure(state='disabled')
									self.logArea.update()

							outfile.write("\n")
							step = step+1
							self.progressbar['value']=int( step*progressBarIncrement)
							self.progressbar.update()



						now = datetime.datetime.now()
						logFile.write("Genotyping sample "+fileRoot1+"  ended at "+now.strftime("%H:%M")+"\n")
						logFile.close()
						outfile.close()
						os.system("rm -f "+dedupFile1+" "+dedupFile2+" *.jf")
						os.system("rm -f alignmenthsbfy43223_sorted.bam null counts.txt")
						os.system("mv *IDCard.txt "+self.OutputFolderEntry.get())
						os.system("mv *_logFile.txt "+self.OutputFolderEntry.get())
				else:
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Please select a valid output folder to start genotyping....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

			else:
				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Please select a valid input file to start genotyping....\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()

		#**************************************************************************
		#***************** Main Genotyping   Algorithm End * **********************
		#**************************************************************************



		#**************************************************************************
		#***************** Main Plotting Algorithm Start  *************************
		#**************************************************************************
		def plotGenotypes():
			genes = ["ul120","ul9","rl12","rl13","ul1","ul11","ul139","rl5a","rl6","ul20","ul146","ul73","ul74"]
			colorDict = {"G1":"blue","G1A":"lightblue","G1B":"darkblue","G1C":"deepskyblue","G2":"red","G2A":"lightcoral","G2B":"darkred","G3":"salmon","G3A":"green","G3B":"darkgreen","G4":"lightgreen","G4A":"yellow","G4B":"darkseagreen","G4C":"pink","G4D":"maroon","G5":"magenta","G6":"orchid","G7":"purple","G8":"silver","G9":"blueviolet","G10":"darkcyan","G11":"navy","G12":"black","G13":"gray","G14":"lightgray"}
			
			if not self.OutputFolderEntry.get() == "Please select a folder....":
				os.system("mkdir "+self.OutputFolderEntry.get())
				packet = StringIO()
				numberOfPlots = tkFileDialog.askopenfilenames(initialdir = "./",title = "Select an input file")
				outputFileName = simpledialog.askstring(title="Output file name",prompt="Insert the output file name")
				

				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Preparing plot....\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()

				shift = 0
				fig = plt.figure(figsize=(10,10))
				plot1 = fig.add_subplot(111)

				for sample in numberOfPlots:
					inputFile = sample
					hg = {}
					infile = open(inputFile)
					while True:
						line = infile.readline().rstrip()
						if not line:
							break
						fields = line.split("\t")
						if not fields[0] in hg:
							hg[fields[0]] = []
						for a in range(1,len(fields),+4):
							hg[fields[0]].append((fields[a],fields[a+1]))

					sizes = []
					colors = []

					for gene in genes:
						for item in hg[gene]:
							colors.append(colorDict[item[0]])
							sizes.append(float(item[1]))
						colors.append("white")
						sizes.append(0.05)


					plot1.pie(sizes,radius=1-shift,colors=colors)
					plot1.pie([1],radius=0.95-shift,colors=["white"])
					
					shift = shift + 0.07

				yLabel = 0.5
				for sample in numberOfPlots:
					plt.text(-0.5,yLabel,(sample.split("/"))[-1],fontsize=14,fontname="arial")
					yLabel = yLabel - 0.07

				legend_elements = []
				for item in colorDict:
					legend_elements.append(Patch(facecolor=colorcode.to_rgba(colorDict[item]),label=item))

				plt.legend(handles=legend_elements,loc=0,bbox_to_anchor=(1.12, 0.9))



				plot1.text(1,0.16,"UL120",rotation=14,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(0.80,0.61,"UL9",rotation=37,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(0.4,0.93,"RL12",rotation=62,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-0.15,1.02,"RL13",rotation=-85,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-0.66,0.84,"UL1",rotation=-60,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-1.05,0.46,"UL11",rotation=-30,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-1.22,-0.05,"UL139",rotation=3,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-1.05,-0.57,"RL5A",rotation=24,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-0.67,-0.94,"RL6",rotation=50,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(-0.17,-1.16,"UL20",rotation=75,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(0.35,-1.13,"UL146",rotation=-70,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(0.75,-0.81,"UL73",rotation=-45,fontweight='bold',fontsize=12,fontname="arial")
				plot1.text(0.97,-0.36,"UL74",rotation=-14,fontweight='bold',fontsize=12,fontname="arial")

				fig.savefig(self.OutputFolderEntry.get()+"/"+outputFileName,dpi=300)

				







				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Plot ready!\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()

				can = None
				new_pdf = None
				existing_pdf = None
				output = None
				page = None
				outputStream = None
			else:
				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Please select a valid output folder to start genotyping....\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()








		def exitProgram():
			exit()







		#'''This class configures and populates the toplevel window.
		#   top is the toplevel containing window.'''
		#_bgcolor = '#d9d9d9'  # X11 color: 'gray85'
		#_fgcolor = '#000000'  # X11 color: 'black'
		#_compcolor = '#d9d9d9' # X11 color: 'gray85'
		#_ana1color = '#d9d9d9' # X11 color: 'gray85'
		#_ana2color = '#ececec' # Closest X11 color: 'gray92'
		#self.style = ttk.Style()
		#if sys.platform == "win32":
	#		self.style.theme_use('winnative')
		#self.style.configure('.',background=_bgcolor)
		#self.style.configure('.',foreground=_fgcolor)
		#self.style.configure('.',font="TkDefaultFont")
		#self.style.map('.',background=
		#	[('selected', _compcolor), ('active',_ana2color)])

		w=840
		h=470
		ws = root.winfo_screenwidth() # width of the screen
		hs = root.winfo_screenheight() # height of the screen
		# calculate x and y coordinates for the Tk root window
		x = (ws/2) - (w/2)
		y = (hs/2) - (h/2)
		top.geometry('%dx%d+%d+%d' % (w, h, x, y))
		top.title("Genotyping")
		top.configure(highlightcolor="black",background="white")

		self.InputFileLabel = ttk.Label(top,background="white",borderwidth=0,foreground="#204949",font="arial")
		self.InputFileLabel.place(x=20, y=20,height=20, width=80)
		self.InputFileLabel.configure(text="Input file")

		self.InputFileEntry = tk.Entry(top,font=("arial",10,"bold"))
		self.InputFileEntry.place(x=20, y=40,height=30, width=300)
		self.InputFileEntry.insert(0,"Please select file....")

		self.InputFileButton = tk.Button(top,command=openInputFile,background="#204949",foreground="white")
		self.InputFileButton.place(x=340,y=40,height=30, width=100)
		self.InputFileButton.configure(text="Open file")


		self.OutputFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949",font=("arial",10,"bold"))
		self.OutputFileLabel.place(x=20, y=100,height=20, width=90)
		self.OutputFileLabel.configure(text="Output folder")

		self.OutputFolderEntry = tk.Entry(top,font=("arial",10,"bold"))
		self.OutputFolderEntry.place(x=20, y=120,height=30, width=300)
		self.OutputFolderEntry.insert(0,"Please select a folder....")

		self.OutputFileButton = tk.Button(top,command=openOutputFolder,background="#204949",foreground="white",font=("arial",10,"bold"))
		self.OutputFileButton.place(x=340,y=120,height=30, width=100)
		self.OutputFileButton.configure(text="Open folder")




		self.dblabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949",font=("arial",10,"bold"))
		self.dblabel.place(x=500, y=20,height=20, width=100)
		self.dblabel.configure(text="kmer database")

		self.dbEntry = tk.Entry(top,font=("arial",10,"bold"))
		self.dbEntry.place(x=500, y=40,height=30, width=200)
		self.dbEntry.insert(0,"mainDB_seqs_filtered.txt")

		self.changeDB = tk.Button(top,background="#204949",foreground="white",font=("arial",10,"bold"))
		self.changeDB.place(x=720, y=40,height=30, width=100)
		self.changeDB.configure(text="Change DB")

		self.CutoffLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949",font=("arial",10,"bold"))
		self.CutoffLabel.place(x=500, y=100,height=20, width=130)
		self.CutoffLabel.configure(text="Detection treshold")

		self.CutoffText = tk.Entry(top, justify='right',font=("arial",10,"bold"))
		self.CutoffText.place(x=500, y=120,height=30, width=150)
		self.CutoffText.insert(0,"0.02")

		self.numThreadsLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949",font=("arial",10,"bold"))
		self.numThreadsLabel.place(x=670,y=100,height=20, width=90)
		self.numThreadsLabel.configure(text="Num. threads")

		self.numThreadsEntry = tk.Entry(top, justify='right',font=("arial",10,"bold"))
		self.numThreadsEntry.place(x=670, y=120,height=30, width=150)
		self.numThreadsEntry.insert(0,"8")

		self.logFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949",font=("arial",10,"bold"))
		self.logFileLabel.place(x=20,y=180,height=20,width=80)
		self.logFileLabel.configure(text="Log window")

		self.logFrame = tk.Frame(top)
		self.logFrame.place(x=20, y=200, height=210, width=600)
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(borderwidth="2")
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(width=125)
		self.logArea = tk.Text(top,state='disabled')
		self.logArea.place(x=25,y=205,height=200, width=590)
		self.logArea.configure(background="white",borderwidth=5)
		self.logArea.configure(selectbackground="#c4c4c4")

		self.progressbar=ttk.Progressbar(top,orient="horizontal",length=600,mode="determinate")
		self.progressbar.place(x=20,y=429)
		self.progressbar['maximum'] = 100

		

		self.runButton = tk.Button(top,command=runGenotyping,background="#204949",foreground="white")
		self.runButton.place(x=650,y=420,height=30,width=80)
		self.runButton.configure(text="Genotype!")

		self.plotButton = tk.Button(top,command=plotGenotypes,background="#204949",foreground="white")
		self.plotButton.place(x=740,y=420,height=30,width=80)
		self.plotButton.configure(text="Plot")

		#self.exitButton = tk.Button(top,command=exitProgram)
		#self.exitButton.place(x=740,y=200,height=30,width=80)
		#self.exitButton.configure(text="Exit")

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Classification.jpg")
		photo1 =ImageTk.PhotoImage(image)
		logoLabel = tk.Label(top, compound=tk.TOP,height=200,width=200,image=photo1,borderwidth=0,highlightthickness=0)
		logoLabel.place(x=640,y=200)
		logoLabel.image = photo1






if __name__ == '__main__':
	vp_start_gui()
