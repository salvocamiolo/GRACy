import tkinter as tk
from tkinter import ttk

from PIL import ImageTk, Image
import sys,os
import tkinter
from tkinter import filedialog as tkFileDialog

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


		def openReferenceFile():
			referenceFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select reference file")
			self.referenceFileEntry.delete(0,tk.END)
			self.referenceFileEntry.insert(0,referenceFile)
			self.referenceFileEntry.icursor(100)

		def openCDSFile():
			cdsFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select CDS file")
			self.cdsFileEntry.delete(0,tk.END)
			self.cdsFileEntry.insert(0,cdsFile)

		def openGFFFile():
			gffFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select GFF file")
			self.gffFileEntry.delete(0,tk.END)
			self.gffFileEntry.insert(0,gffFile)

		def openPos2PlotFile():
			pos2plotFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
			self.posToPlotEntry.delete(0,tk.END)
			self.posToPlotEntry.insert(0,pos2plotFile)

		def openInputFile():
			inputFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select input file")
			self.inputFileEntry.delete(0,tk.END)
			self.inputFileEntry.insert(0,inputFile)

		def openOutputFolder():
			outputolder = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
			self.outputFolderEntry.delete(0,tk.END)
			self.outputFolderEntry.insert(0,outputolder)


		def mainAlgorithm():
			if self.referenceFileEntry.get() == "Please select a file...." or self.inputFileEntry.get()=="Please select a file...." or self.gffFileEntry.get() =="Please select a file...." or self.cdsFileEntry.get()=="Please select a file...." or self.outputFolderEntry.get()=="Please select a folder....":
				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "Please enter a valid input, reference, cds, annotation files and output folder  to continue!\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()

			else:
				infile = open(self.inputFileEntry.get())
				line = infile.readline().rstrip()
				suffix = (line.split("."))[-1]
				numThreads = self.numThreadsEntry.get()
				outputFolder = self.outputFolderEntry.get()
				referenceFile = self.referenceFileEntry.get()
				if suffix == "fastq" or suffix == "fq":
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "The input file contains fastq reads!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()


					#***********************************************************************
					#********************* Main algorithm reads alignment start ************
					#***********************************************************************

					infile.close()

					infile = open(self.inputFileEntry.get())

					while True:
						read1 = infile.readline().rstrip()
						if not read1:
							break

						read2 = infile.readline().rstrip()
						if not read1:
							break

						if ".fastq" in read1:
							sampleName = ((read1.split("/"))[-1].split("_1.fastq"))[0]
						if ".fq" in read1:
							sampleName = ((read1.split("/"))[-1].split("_1.fq"))[0]

						sampleTempFolder = sampleName+"_4864985345_tempFolder"
						os.system("mkdir -p "+sampleTempFolder)
						os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build "+referenceFile+" ./"+sampleTempFolder+"/reference -q")

						print("Analising sample",sampleName)
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "Analising sample"+sampleName+"....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Reads quality filtering before alignment\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Running Trimgalore....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()



						print("Reads quality filtering before alignment")
						print("Running Trimgalore")
						prefix1 = ((read1.split("/"))[-1]).replace(".fastq","")
						prefix2 = ((read2.split("/"))[-1]).replace(".fastq","")
						os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt --paired -q 0  "+read1+" "+read2 + " -o "+sampleTempFolder)

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Performing deduplication....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						print("Performing deduplication")
						os.system("echo "+sampleTempFolder+"/"+prefix1+"_val_1.fq > "+sampleTempFolder+"/inputFastUniq")
						os.system("echo "+sampleTempFolder+"/"+prefix2+"_val_2.fq >> "+sampleTempFolder+"/inputFastUniq")
						os.system(installationDirectory+"src/conda/bin/fastuniq -i "+sampleTempFolder+"/inputFastUniq -t q -o "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_1.fastq -p "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_2.fastq")

						print("Performing prinseq quality filtering")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Running Prinseq....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/prinseq-lite.pl -fastq "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_1.fastq  -fastq2 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_2.fastq -min_qual_mean 25 -trim_qual_right 30 -trim_ns_right 20  -trim_qual_window 5 -trim_qual_step 1 -min_len 80 -out_bad null -out_good "+sampleName+"_trimmed_dedup_pr")
						os.system("mv "+sampleName+"_trimmed_dedup_pr* ./"+sampleTempFolder)
						os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/trimPolyN.py "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_1.fastq "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_2.fastq")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Reads alignment on reference\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Performing alignment....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						print("Aigning reads to reference")
						os.system(installationDirectory+"src/conda/bin/bowtie2 --end-to-end -1 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_1.fastq -2 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_2.fastq -x "+sampleTempFolder+"/reference -S "+sampleTempFolder+"/"+sampleName+"_alignment.sam -p "+numThreads)

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Converting sam to bam....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						print("Converting sam to bam")
						os.system(installationDirectory+"src/conda/bin/samtools view -bS -h "+sampleTempFolder+"/"+sampleName+"_alignment.sam > "+sampleTempFolder+"/"+sampleName+"_alignment.bam")

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Sorting bam....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						print("Sorting bam")
						os.system(installationDirectory+"src/conda/bin/samtools sort -o "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+sampleTempFolder+"/"+sampleName+"_alignment.bam")

						if self.dedupChkValue.get()==True:
							self.logArea.configure(state='normal')
							self.logArea.insert(tk.END, "*  *  Marking duplicates....\n")
							self.logArea.see(tk.END)
							self.logArea.configure(state='disabled')
							self.logArea.update()

							print("Marking duplicates ")
							os.system(installationDirectory+"src/conda/bin/picard  AddOrReplaceReadGroups I="+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam O="+sampleTempFolder+"/"+sampleName+"_rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")
							os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I="+sampleTempFolder+"/"+sampleName+"_rg_added_sorted.bam O="+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M="+sampleTempFolder+"/output.metrics")

						else:
							os.system("mv "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam")
						print("Calling snps with lofreq")
						#Analyze snps with lowfreq
						os.system(installationDirectory+"/src/conda/bin/samtools faidx "+referenceFile)
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Calling SNPs with lofreq\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()

						os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.numThreadsEntry.get()+" -q "+str(self.BaseQualEntry.get())+" -Q "+str(self.BaseQualEntry.get()) +" -f "+referenceFile+" -o "+sampleTempFolder+"/"+sampleName+"_SNPs.vcf "+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam")
						print("Lofreq terminated")
						#sys.stdin.read(1)

						#Analyze indes using the GATK pipeline

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Calling INDELs with GATK\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()


						os.system(installationDirectory+"src/conda/bin/samtools faidx "+referenceFile)

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  *  Calling INDELs....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						#os.system("java -jar "+installationDirectory+"resources/picard.jar CreateSequenceDictionary R= "pfal.fa O= pfal.fa.dict")
						#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R "+referenceFile+" -I "+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam  -o "+sampleTempFolder+"/"+sampleName+"_output.vcf -A StrandAlleleCountsBySample")
						#os.system("mv "+sampleTempFolder+"/"+sampleName+"_output.vcf "+sampleTempFolder+"/"+sampleName+"_indels.vcf")


						os.system("mv "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+outputFolder+"/")
						os.system("rm "+sampleTempFolder+"/"+sampleName+"_alignment* -f")
						#os.system("rm -f inputFastUniq "+prefix1+"_val_1.fq "+prefix2+"_val_2.fq  "+sampleName+"_trimmed_dedup_* "+sampleName+"_markedDuplicates.ba* "+sampleName+"_rg_added_sorted.bam "+sampleName+"*vcf.idx*")
						os.system(" ls "+sampleTempFolder+"/*SNPs.vcf > "+sampleTempFolder+"/vcfFilesToExamine")
						os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/polyAn.py "+sampleTempFolder+"/vcfFilesToExamine"+" "+self.cdsFileEntry.get()+" "+self.gffFileEntry.get())
						#os.system("rm vcfFilesToExamine *.dict *.fai *trimming_report.txt")
						os.system("ls "+sampleTempFolder+"/*_snpEffect.txt > "+sampleTempFolder+"/file2Plot")
						os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/buildTable.py "+sampleTempFolder+"/file2Plot")
						os.system("mv  "+sampleTempFolder+"/*_snpEffect.txt "+sampleTempFolder+"/*_snpFreq.txt "+sampleTempFolder+"/snpTable.txt "+sampleTempFolder+"/*.vcf "+outputFolder+"/" )
						os.system("mv "+sampleTempFolder+"/* ./"+outputFolder+"/")
						os.system("rm -rf "+sampleTempFolder)


				#Not used
				if suffix == "vcf":
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "The input file contains vcf datasets!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/polyAn.py "+self.inputFileEntry.get()+" "+self.cdsFileEntry.get()+" "+self.gffFileEntry.get())
					os.system("rm vcfFilesToExamine *.dict  mapped.bam  rg_added_sorted.bam *trimming_report.txt")
					os.system("ls *_snpEffect.txt > file2Plot")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/buildTable.py file2Plot")
					os.system("mv  *_snpEffect.txt *_snpFreq.txt snpTable.txt "+outputFolder+"/" )




				if not suffix == "fastq" and not suffix == "vcf" and not suffix == "fq":
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "The input file does not contain datasets GRACy can recognise!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()



		def exitProgram():
			exit()

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
		h=520
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
		self.inputFileLabel.place(x=20,y=20,width=60,height=20)

		self.inputFileEntry = tk.Entry(top)
		self.inputFileEntry.place(x=20,y=40,width=750,height=30)
		self.inputFileEntry.insert(0,"Please select a file....")

		self.inputFileButton = tk.Button(top,command=openInputFile,background="#204949",foreground="white")
		self.inputFileButton.place(x=780,y=40,width=100,height=30)
		self.inputFileButton.configure(text="Open file")

		#self.posToPlotLabel = tk.Label(top)
		#self.posToPlotLabel.configure(text="Positions to plot")
		#self.posToPlotLabel.place(x=470,y=20,width=110,height=20)

		#self.posToPlotEntry = tk.Entry(top)
		#self.posToPlotEntry.place(x=470,y=40,width=300,height=30)
		#self.posToPlotEntry.insert(0,"Please select a file....")

		#self.posToPlotButton = tk.Button(top,command=openPos2PlotFile)
		#self.posToPlotButton.place(x=780,y=40,width=100,height=30)
		#self.posToPlotButton.configure(text="Open file")


		self.cdsFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.cdsFileLabel.configure(text="CDS file")
		self.cdsFileLabel.place(x=20,y=80,width=50,height=20)

		self.cdsFileEntry = tk.Entry(top)
		self.cdsFileEntry.place(x=20,y=100,width=300,height=30)
		self.cdsFileEntry.insert(0,"Please select a file....")

		self.cdsFileButton = tk.Button(top,command=openCDSFile,background="#204949",foreground="white")
		self.cdsFileButton.place(x=330,y=100,width=100,height=30)
		self.cdsFileButton.configure(text="Open file")


		self.gffFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.gffFileLabel.configure(text="Annotation file")
		self.gffFileLabel.place(x=470,y=80,width=100,height=20)

		self.gffFileEntry = tk.Entry(top)
		self.gffFileEntry.place(x=470,y=100,width=300,height=30)
		self.gffFileEntry.insert(0,"Please select a file....")

		self.gffFileButton = tk.Button(top,command=openGFFFile,background="#204949",foreground="white")
		self.gffFileButton.place(x=780,y=100,width=100,height=30)
		self.gffFileButton.configure(text="Open file")


		self.referenceFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.referenceFileLabel.configure(text="Reference file")
		self.referenceFileLabel.place(x=20,y=140,width=95,height=20)

		self.referenceFileEntry = tk.Entry(top)
		self.referenceFileEntry.place(x=20,y=160,width=300,height=30)
		self.referenceFileEntry.insert(0,"Please select a file....")

		self.referenceFileButton = tk.Button(top,command=openReferenceFile,background="#204949",foreground="white")
		self.referenceFileButton.place(x=330,y=160,width=100,height=30)
		self.referenceFileButton.configure(text="Open file")


		self.outputFolderLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.outputFolderLabel.configure(text="Output folder")
		self.outputFolderLabel.place(x=470,y=140,width=95,height=20)

		self.outputFolderEntry = tk.Entry(top)
		self.outputFolderEntry.place(x=470,y=160,width=300,height=30)
		self.outputFolderEntry.insert(0,"Please select a folder....")

		self.outputFolderButton = tk.Button(top,command=openOutputFolder,background="#204949",foreground="white")
		self.outputFolderButton.place(x=780,y=160,width=100,height=30)
		self.outputFolderButton.configure(text="Open folder")


		#self.freqCutoffLabel = tk.Label(top)
		#self.freqCutoffLabel.configure(text="Freq cutoff")
		#self.freqCutoffLabel.place(x=20,y=220,width=80,height=20)

		#self.freqCutoffEntry = tk.Entry(top,justify='right')
		#self.freqCutoffEntry.place(x=20,y=240,width=100,height=30)
		#self.freqCutoffEntry.insert(0,"0.1")

		self.numThreadsLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.numThreadsLabel.configure(text="Num. threads")
		self.numThreadsLabel.place(x=20,y=220,width=90,height=20)

		self.numThreadsEntry = tk.Entry(top,justify='right')
		self.numThreadsEntry.place(x=20,y=240,width=100,height=30)
		self.numThreadsEntry.insert(0,"8")

		self.BaseQualLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.BaseQualLabel.configure(text="Base Qual. cutoff")
		self.BaseQualLabel.place(x=150,y=220,width=120,height=20)

		self.BaseQualEntry = tk.Entry(top,justify='right')
		self.BaseQualEntry.place(x=150,y=240,width=120,height=30)
		self.BaseQualEntry.insert(0,"30")



		self.dedupChkValue = tk.BooleanVar()
		self.dedupChkValue.set(True)
		self.dedupCheckButton = tk.Checkbutton(top,variable=self.dedupChkValue,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.dedupCheckButton.place(x=280,y=250,height=20,width=180)
		self.dedupCheckButton.configure(text="Perform deduplication")

		self.runButton = tk.Button(top,command=mainAlgorithm,background="#204949",foreground="white")
		self.runButton.place(x=700,y=470,width=180,height=30)
		self.runButton.configure(text="Run")

		#self.plotButton = tk.Button(top)
		#self.plotButton.place(x=660,y=240,width=100,height=30)
		#self.plotButton.configure(text="Plot")

		#self.exitButton = tk.Button(top,command=exitProgram)
		#self.exitButton.place(x=660,y=240,width=100,height=30)
		#self.exitButton.configure(text="Exit")

		self.logFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.logFileLabel.place(x=20,y=280,height=20,width=100)
		self.logFileLabel.configure(text="Log window")

		self.logFrame = tk.Frame(top)
		self.logFrame.place(x=20, y=280, height=220, width=660)
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(borderwidth="2")
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(width=125)
		self.logArea = tk.Text(top,state='disabled')
		self.logArea.place(x=25,y=285,height=210, width=650)
		self.logArea.configure(background="white",borderwidth=5)
		self.logArea.configure(selectbackground="#c4c4c4")

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/SNPcalling.jpg")
		photo1 =ImageTk.PhotoImage(image)
		logoLabel = tk.Label(top, compound=tk.TOP,height=200,width=200,image=photo1,borderwidth=0,highlightthickness=0)
		logoLabel.place(x=690,y=260)
		logoLabel.image = photo1




if __name__ == '__main__':
	vp_start_gui()
