import tkinter as tk
from tkinter import ttk

from PIL import ImageTk, Image

import sys,os
import datetime
import time
from Bio import SeqIO

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
	#rt = root
	w = tk.Toplevel (root)
	#VirHosFilt_support.set_Tk_var()
	top = Toplevel1 (w)

	return (w, top)

def destroy_Toplevel1():
	global w
	w.destroy()
	w = None

class Toplevel1:
	def __init__(self, top=None):

		def bowtiePE(reference,read1,read2,numTh):
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build "+reference+" reference -q >null 2>&1")

			os.system(installationDirectory+"src/conda/bin/bowtie2  -1 "+read1+" -2 "+read2+" -x reference -S test.sam -p "+numTh)
			os.system(installationDirectory+"src/conda/bin/samtools view -bS -h test.sam > test.bam 2>null")
			os.system(installationDirectory+"src/conda/bin/samtools sort -o test_sorted.bam test.bam >null 2>&1")



		def bwaPE(reference,read1,read2,alName,numThreads,editDist):
			os.system(installationDirectory+"src/conda/bin/bwa index "+reference+" >null 2>&1")
			os.system(installationDirectory+"src/conda/bin/bwa aln "+reference+" "+read1+" -t "+numThreads+" -n "+editDist+" -k "+editDist+" >read1.sai 2>null")
			os.system(installationDirectory+"src/conda/bin/bwa aln "+reference+" "+read2+" -t "+numThreads+" -n "+editDist+" -k "+editDist+" >read2.sai 2>null")
			os.system(installationDirectory+"src/conda/bin/bwa sampe "+reference+" read1.sai read2.sai "+read1+" "+read2+" >"+alName+".sam 2>null")
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/cleanSoftAndUnmapped.py "+alName+".sam")
			os.system(installationDirectory+"src/conda/bin/samtools view -bS "+alName+".sam_cleaned.sam > "+alName+".bam 2>null")
			os.system(installationDirectory+"src/conda/bin/samtools sort -o "+alName+"_sorted.bam "+alName+".bam >null 2>&1")
			os.system(installationDirectory+"src/conda/bin/samtools index "+alName+"_sorted.bam >null 2>&1")

		def exitProgram():
			exit()

		def openInputFile():
			inputFile = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
			self.inputFileEntry.delete(0,tk.END)
			self.inputFileEntry.insert(0,inputFile)

		#***************************************************
		#*************** Start main algorithm **************
		#***************************************************

		def performAssembly():
			infileName = self.inputFileEntry.get()

			infile = open(infileName)
			while True:
				line = infile.readline().rstrip()
				if not line:
					break
				confFile = open(line)
				projectName = ((confFile.readline().rstrip()).split("\t"))[1]
				os.system("mkdir -p "+projectName)
				os.system("cp "+line+" "+projectName)
				os.chdir(projectName)
				logFile = open('logFile.log','w')
				read1 = ((confFile.readline().rstrip()).split("\t"))[1]
				read2 = ((confFile.readline().rstrip()).split("\t"))[1]
				read1_toFill = ((confFile.readline().rstrip()).split("\t"))[1]
				read2_toFill = ((confFile.readline().rstrip()).split("\t"))[1]
				confFile.readline() #Read comment

				now = datetime.datetime.now()
				logFile.write("Date: "+now.strftime("%Y-%m-%d")+"\n")
				logFile.write("Sample name: "+projectName+"\n")
				logFile.write("Read1 fastq: "+read1+"\n")
				logFile.write("Read2 fastq: "+read2+"\n\n\n")

				#***************************************************************
				#**************** 1 Reads quality filtering ********************
				#***************************************************************


				qualityFiltering = ((confFile.readline().rstrip()).split("\t"))[1]
				logFile.write("Quality filtering: "+qualityFiltering+"\n")

				if qualityFiltering == "yes" or qualityFiltering == "Yes":
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Starting filtering for sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					now = datetime.datetime.now()
					logFile.write("Quality filtering started at "+now.strftime("%H:%M")+"\n")

					qualConfFile = open("qualityFiltering.conf","w")
					qualConfFile.write("ProjectName\tallSamples\n")
					qualConfFile.write("Sample_start****************************************"+"\n")
					qualConfFile.write("SampleName\t"+projectName+"\n")
					qualConfFile.write("Read1\t"+read1+"\n")
					qualConfFile.write("Read2\t"+read2+"\n")
					qualConfFile.write("Sample_start***Filtering and trimming options\n")
					for a in range(6):
						qualConfFile.write(confFile.readline())
					qualConfFile.write("SampleEnd*******************************************\n")
					qualConfFile.close()
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Running prinseq.....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/runQualityFiltering.py qualityFiltering.conf "+installationDirectory)
					os.system("mkdir -p 1_cleanReads")
					os.system("mv "+projectName+"_hq_1.fastq ./1_cleanReads/qualityFiltered_1.fq")
					os.system("mv "+projectName+"_hq_2.fastq ./1_cleanReads/qualityFiltered_2.fq")


					now = datetime.datetime.now()
					logFile.write("Quality filtering ended at "+now.strftime("%H:%M")+"\n\n")
					os.system("rm -f badReads1.fastq badReads2.fastq *filterStats.txt *singletons.fastq paired_normalized.fq paired.fq qualityFiltering.conf")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Filtering completed for sample "+projectName+"!!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

				else:
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "No quality filtering required on project  "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					for a in range(6):
						confFile.readline()
					if os.path.isdir("./1_cleanReads/") == True:
						if os.path.isfile("./1_cleanReads/qualityFiltered_1.fq") == False:
							print("File ","./1_cleanReads/qualityFiltered_1.fq does not exist, now exiting....")
							exit()
						if os.path.isfile("./1_cleanReads/qualityFiltered_2.fq") == False:
							print("File ","./1_cleanReads/qualityFiltered_2.fq does not exist, now exiting....")
							exit()
					else:
						os.system("mkdir 1_cleanReads")
						os.chdir("1_cleanReads")
						os.system("ln -s "+read1+" qualityFiltered_1.fq")
						os.system("ln -s "+read2+" qualityFiltered_2.fq")
						os.chdir("../")


				#***************************************************************
				#*********************** 2 Denovo assembly *********************
				#***************************************************************
				confFile.readline() #Read comment

				denovoAssembly = ((confFile.readline().rstrip()).split("\t"))[1]
				print(denovoAssembly)
				if denovoAssembly == "yes" or denovoAssembly=="Yes":
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nPerforming denovo assembly on sample "+projectName+"\n")
					self.logArea.insert(tk.END, "*  Normalizing reads....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					now = datetime.datetime.now()
					logFile.write("De novo assembly started at "+now.strftime("%H:%M")+"\n")
					print("\nPerforming denovo assembly.......")
					os.system(installationDirectory+"src/conda/bin/interleave-reads.py 1_cleanReads/qualityFiltered_1.fq 1_cleanReads/qualityFiltered_2.fq -o paired.fq")
					availableMemory = self.memoryEntry.get()
					os.system(installationDirectory+"src/conda/bin/normalize-by-median.py  -k 17 -C 200 -M "+availableMemory+"e9 -p  -o - paired.fq > paired_normalized.fq")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/splitIntervealed.py paired_normalized.fq")
					os.system("mv newRead_1.fastq ./1_cleanReads/"+projectName+"_hq_1.fastq")
					os.system("mv newRead_2.fastq ./1_cleanReads/"+projectName+"_hq_2.fastq")
					os.system("rm -f paired_normalized.fq paired.fq")

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Running SPAdes....")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					now = datetime.datetime.now()
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getBestAssembly.py ./1_cleanReads/"+projectName+"_hq_1.fastq ./1_cleanReads/"+projectName+"_hq_2.fastq assemblyStatistics.txt "+installationDirectory)
					os.system("mkdir 2_spadesAssembly")
					os.system("mv scaffolds.fasta ./2_spadesAssembly/")
					os.system("rm numReads.txt subsample_1.fq subsample_2.fq null N50.txt")
					os.system("mv assemblyStatistics.txt ./2_spadesAssembly")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "De novo assembly completed on sample "+projectName+"!!\n")
					self.logArea.insert(tk.END, "Normalizing reads....")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()


				else:
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nNo de novo required for sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					if os.path.isfile("./2_spadesAssembly/scaffolds.fasta") == False:
						print("You chose not to run assembler but scaffolds.fasta file is not there. Now exiting......")
						logFile.write("ERROR!\n You chose not to run assembler but scaffolds.fasta file is not there\n")
						exit()
				if os.path.isfile("./2_spadesAssembly/scaffolds.fasta") == False:
					print("Something went wrong with the assembly. Now exiting......")
					logFile.write("ERROR!\n Something went wrong with the assembly\n")
					exit()
				else:
					os.system("rm -rf ./2_spadesAssembly/corrected")
					now = datetime.datetime.now()
					logFile.write("De novo assembly ended at "+now.strftime("%H:%M")+"\n\n")



				#***************************************************************
				#******************* 3 Scaffold Oriantation ********************
				#***************************************************************
				confFile.readline() #Read comment
				performScaffolding = ((confFile.readline().rstrip()).split("\t"))[1]

				if performScaffolding == "yes" or performScaffolding == "Yes":
					now = datetime.datetime.now()
					logFile.write("Scaffolding started at "+now.strftime("%H:%M")+"\n\n")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nPerforming scaffodling on sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("mkdir -p 3_scaffoldsOrientation")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/retrieveNodes.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/getSequenceFromFasta.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds_careful.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds_trivial.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/utils/biomodule.py ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/sequence_a.fasta ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/sequence_a_prime.fasta ./3_scaffoldsOrientation/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/createCenterScaffold.py ./3_scaffoldsOrientation/")
					os.system("cp "+installationDirectory+"data/merlinReference/merlinGenome_190000_200000_f.txt ./3_scaffoldsOrientation/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py ./3_scaffoldsOrientation/")
					os.system("cp "+installationDirectory+"/src/scripts/assembly/utils/gapPrediction.py ./3_scaffoldsOrientation/")
					os.system("cp "+installationDirectory+"data/merlinReference/hcmv_genomes.fasta* ./3_scaffoldsOrientation/")
					#os.system("cp "+installationDirectory+"src/scripts/assembly/utils/cap3 ./3_scaffoldsOrientation/")

					os.chdir("3_scaffoldsOrientation")
					print("\nCalculating the average insert size......")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Calculating the average insert size\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					bowtiePE("../2_spadesAssembly/scaffolds.fasta","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())
					os.system(installationDirectory+"src/conda/bin/picard CollectInsertSizeMetrics I=test_sorted.bam  O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5")
					os.system("head -8 insert_size_metrics.txt | tail -2 | cut -f 6 | tail -1 >insert.size")
					isize = open("insert.size")
					insertSize = isize.readline().rstrip()
					isize.close()
					print("insertSize",insertSize)
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Insert size = "+insertSize+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					print("\nAligning the contigs to the merlin reference gneome......")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Aligning contigs to reference\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()


					os.system(installationDirectory+"src/conda/bin/python createCenterScaffold.py "+projectName+" "+installationDirectory+" "+insertSize+" "+self.threadsEntry.get())
					os.system("cp newFinalScaffold.fasta longestScaffold.fasta ")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Filling gaps from existing sequences\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python gapPrediction.py longestScaffold.fasta ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq "+installationDirectory)

					finalScaffoldFile = open("finalScaffold.fasta","w")
					finalScaffoldFile.write(">finalScaffold\n")
					gapPresent = 0
					for seq_record in SeqIO.parse("filledGenome.fasta","fasta"):
						finalScaffoldFile.write(str(seq_record.seq)+"\n")
						if "N" in str(seq_record.seq):
							gapPresent = 1
					finalScaffoldFile.close()
					if gapPresent == 1:
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Filling gaps with existing sequences\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						gfFile = open("gapfillerlib.txt","w")
						gfFile.write("lib1 bwa ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq "+ insertSize+" 0.25 FR")
						gfFile.close()
						os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/GapFiller -q "+installationDirectory+" -l gapfillerlib.txt -s filledGenome.fasta -T "+self.threadsEntry.get())
						os.system("mv finalScaffold.fasta finalScaffold.fasta_beforeGapfilling")
						os.system("cp standard_output/standard_output.gapfilled.final.fa ./finalScaffold.fasta")

					now = datetime.datetime.now()
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Scaffolding completed on sample "+projectName+"!!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					logFile.write("Scaffolding finishes at "+now.strftime("%H:%M")+"\n")
					if self.verboseChkValue.get() == False:
						os.system("rm -rf *.bam *.sam *.bt2 *txt *ap* cent* filled* hcmv* join* last* long* merlin* null *fastq sb* sequence* standard* temp* two* *.py found* new*")
					os.chdir("../")
				else:
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nNo scaffolding required for sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					if os.path.isfile("./3_scaffoldsOrientation/finalScaffold.fasta") == False:
						print("You chose not to run the scaffolding step but finalScaffold.fasta file is not there. Now exiting......")
						logFile.write("You chose not to run the scaffolding step but finalScaffold.fasta file is not there. Now exiting......")
						exit()



				#***************************************************************
				#********************** 4 Create Consensus *********************
				#***************************************************************
				confFile.readline() #Read comment
				perform1stConsensusCalling = ((confFile.readline().rstrip()).split("\t"))[1]
				if perform1stConsensusCalling == "yes" or perform1stConsensusCalling == "Yes":
					now = datetime.datetime.now()
					logFile.write("First consensus calling started at "+now.strftime("%H:%M")+"\n\n")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nFirst round of consensus calling on sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("mkdir 4_createConsensus")
					os.system("cp ./3_scaffoldsOrientation/finalScaffold.fasta ./4_createConsensus")
					#os.system("cp "+installationDirectory+"assembly/scripts/joinConsensus.py ./4_createConsensus")
					#os.system("cp "+installationDirectory+"assembly/scripts/hcmvConsensusCallPipeline ./4_createConsensus")
					#os.system("cp "+installationDirectory+"assembly/scripts/createConsensus ./4_createConsensus")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py ./4_createConsensus")
					#os.system("cp "+installationDirectory+"assembly/scripts/libdeflate.so ./4_createConsensus")
					os.chdir("4_createConsensus")

					for seq_record in SeqIO.parse("finalScaffold.fasta","fasta"):
						assemblyLength = len(str(seq_record.seq))


					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing first portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python ./extractSeqByRange.py finalScaffold.fasta finalScaffold 1 15001 f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					bowtiePE("finalScaffold_1_15001_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()


					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")# >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_1_15001_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_1_15001_f.txt")

					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_1_15001_f.txt -o output.vcf dedupped.bam")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_1_15001_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_1_15001_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_1_15001_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_firstPortion.vcf")
					os.system("mv output_filtered.vcf  output_firstPortion_filtered.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")
					


					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing second portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python ./extractSeqByRange.py finalScaffold.fasta finalScaffold 15001 "+str(assemblyLength -10000 )+" f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					bowtiePE("finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt")
					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt -o output.vcf dedupped.bam")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					#os.system(installationDirectory+"resources/bgzip -c output.vcf_filtered.vcf > output.vcf_filtered.vcf.gz 2>null")
					#os.system(installationDirectory+"resources/tabix output.vcf_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_secondPortion.vcf")
					os.system("mv output_filtered.vcf  output_filtered_secondPortion.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")


					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing third portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python ./extractSeqByRange.py finalScaffold.fasta finalScaffold "+str(assemblyLength - 10000 )+" 2000000 f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					bowtiePE("finalScaffold_"+str(assemblyLength -10000 )+"_2000000_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt")
					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt -o output.vcf dedupped.bam")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils//vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					#os.system(installationDirectory+"resources/filterVCF.py output.vcf >null 2>&1")
					#os.system(installationDirectory+"resources/bgzip -c output.vcf_filtered.vcf > output.vcf_filtered.vcf.gz 2>null")
					#os.system(installationDirectory+"resources/tabix output.vcf_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_thirdPortion.vcf")
					os.system("mv output_filtered.vcf  output_filtered_thirdPortion.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")



					finalSequence = ""

					for seq_record in SeqIO.parse("finalScaffold_1_15001_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					for seq_record in SeqIO.parse("finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					for seq_record in SeqIO.parse("finalScaffold_"+str(assemblyLength -10000 )+"_2000000_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					outfile = open(projectName+"_genome.fasta","w")
					outfile.write(">finalScaffold\n"+finalSequence+"\n")




					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "First round of consensus calling completed on sample "+projectName+"!!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					now = datetime.datetime.now()
					logFile.write("First consensus calline ended at "+now.strftime("%H:%M")+"\n\n")
					os.chdir("../")
				else:
					if os.path.isfile("./4_createConsensus/"+projectName+"_genome.fasta")==False:
						print("You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......")
						logFile.write("You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......")
						exit()



				#***************************************************************
				#********************** 5 Refine Consensus *********************
				#***************************************************************
				confFile.readline() #Read comment
				refineAssembly = ((confFile.readline().rstrip()).split("\t"))[1]
				if refineAssembly == "yes" or refineAssembly == "Yes":
					now = datetime.datetime.now()
					logFile.write("Assembly sequence refining started at "+now.strftime("%H:%M")+"\n\n")

					os.system("mkdir -p 5_refineAssembly")
					os.system("cp ./4_createConsensus/"+projectName+"_genome.fasta ./5_refineAssembly/finalScaffold.fasta")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds_careful.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/joinScaffolds_trivial.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/utils/biomodule.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/maskLowCoverage.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/cleanSoftAndUnmapped.py ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"data/merlinReference/sequence_a.fasta ./5_refineAssembly/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/revComp.py ./5_refineAssembly/")
					os.chdir("5_refineAssembly")


			

					#Start the refine assembly algorithm

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nRefining assembly on sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					notRefined = open("notRefined.txt","w")


					os.system("cp finalScaffold.fasta newReference.fasta")
					for seq_record in SeqIO.parse("finalScaffold.fasta","fasta"):
								reference = str(seq_record.seq)

					newBits = []
					previousEnd = []
					previousEnd.append(0)
					eof = 0
					while True and eof == 0:
						for seq_record in SeqIO.parse("newReference.fasta","fasta"):
								reference = str(seq_record.seq)

						print("Performing Alignment....")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Aligning reads on assembly\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						bwaPE("newReference.fasta",read1,read2,"test",self.threadsEntry.get(),"0.02") #Change here the number of threads
						os.system(installationDirectory+"src/conda/bin/samtools faidx newReference.fasta >null 2>&1")
						print("Calculating coverage on assembly....")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Looking for low coverage regions \n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"src/conda/bin/samtools mpileup -f newReference.fasta test_sorted.bam >finalPileup.txt 2>null")
						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Masking low coverage regions \n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						os.system(installationDirectory+"/src/conda/bin/python ./maskLowCoverage.py newReference.fasta finalPileup.txt ")
						ranges = open("Nranges_smooth.txt")
						line = ranges.readline().rstrip()
						if not line:
							break
						fields = line.split("\t")

						while int(fields[0]) <= previousEnd[len(previousEnd)-1]+50:
							line = ranges.readline().rstrip()
							if not line:
								eof = 1
								break
							fields = line.split("\t")

						previousEnd.append(int(fields[1]))
						print("Refining the range",line)
						print("Length assembly",len(reference),"From",fields[0],"To",fields[1])

						self.logArea.configure(state='normal')
						self.logArea.insert(tk.END, "*  Filling region "+line+"....\n")
						self.logArea.see(tk.END)
						self.logArea.configure(state='disabled')
						self.logArea.update()
						if len(line)>2 and int(fields[0])>3000 and int(fields[0])< (len(reference) - 3000) and  int(fields[1])>3000 and int(fields[1])< (len(reference) - 3000):
							os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[0])-1500)+" "+str(int(fields[0])-500)+" f")
							os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py newReference.fasta finalScaffold "+str(int(fields[1])+500)+" "+str(int(fields[1])+1500)+" f")

							numTries = 0
							while True:
								numTries += 1
								if numTries ==2:
									break
								print("Performing joinScaffold_careful algorithm on range",line,"....")
								os.system(installationDirectory+"src/conda/bin/python joinScaffolds_careful.py join "+read1_toFill+"  "+read2_toFill+"  finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f "+installationDirectory+" "+self.threadsEntry.get())
								if os.path.isfile("joined_finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==False:
									print("Refining range",fields[0],fields[1],"with joinScaffold_careful failed. Now trying joinScaffold on range",line,"....")
									os.system(installationDirectory+"src/conda/bin/python joinScaffolds.py join "+read1_toFill+"  "+read2_toFill+"  finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f "+installationDirectory+" "+self.threadsEntry.get())
									if os.path.isfile("joined_finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==False:
										print("Both alogrithms failed on range,",line)
										print("Now performing 30 steps of joinScaffold trivial in both directions and recording the output")
										os.system(installationDirectory+"src/conda/bin/python joinScaffolds_trivial.py "+ read1_toFill+" "+read2_toFill+" finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt f finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt f 30 "+installationDirectory )
										os.system("mv joinScaffold_trivialSeq.fasta finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_greedy.fasta ")
										os.system(installationDirectory+"src/conda/bin/python joinScaffolds_trivial.py "+ read1_toFill+" "+read2_toFill+" finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt r finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt r 30 "+installationDirectory )
										os.system(installationDirectory+"src/conda/bin/python revComp.py joinScaffold_trivialSeq.fasta >temp.fasta; mv temp.fasta joinScaffold_trivialSeq.fasta")
										os.system("mv joinScaffold_trivialSeq.fasta finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_greedy.fasta")

									else:
										print("The algorithm joinScaffold was sucessful on range",line)
										break
								else:
									print("The algorithm joinScaffold_careful was sucessful on range",line)
									break


							if os.path.isfile("joined_"+"finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt")==True:
								for seq_record in SeqIO.parse("joined_"+"finalScaffold_"+str(int(fields[0])-1500)+"_"+str(int(fields[0])-500)+"_f.txt_finalScaffold_"+str(int(fields[1])+500)+"_"+str(int(fields[1])+1500)+"_f.txt","fasta"):
									newBit = str(seq_record.seq)
								if len(newBit) > 1600:
									newReference = reference[:reference.find(newBit[:100])]+newBit+reference[reference.find(newBit[-100:])+100:]
									newRef = open("newReference.fasta","w")
									newRef.write(">finalScaffold\n"+newReference+"\n")
									newRef.close()


					noRef = []
					infile2 = open("notRefined.txt")
					while True:
						line = infile2.readline()
						if not line:
							break
						noRef.append(line)
					infile2.close()
					if len(noRef) >0:
						logFile.write("WARNING! The following ranges were not closed during the assembly refining step:\n")
						for line in noRef:
							logFile.write(line)

					for seq_record in SeqIO.parse("newReference.fasta","fasta"):
						newRefSeq = str(seq_record.seq)
					finalScaffoldFile = open("finalScaffold.fasta","w")
					finalScaffoldFile.write(">finalScaffold\n"+newRefSeq+"\n")
					finalScaffoldFile.close()


					now = datetime.datetime.now()
					logFile.write("Assembly sequence refining ended at "+now.strftime("%H:%M")+"\n")

					os.chdir("../")

				else:
					if os.path.isfile("./5_refineAssembly/finalScaffold.fasta")==False:
						print("You chose not to run the assembly refining step but the finalScaffold.fasta file is not there. Now exiting......")
						logFile.write("You chose not to run the assembly refining step but the finalScaffold.fasta file is not there. Now exiting......")
						exit()




				#***************************************************************
				#********************** 6 Create Consensus *********************
				#***************************************************************
				confFile.readline() #Read comment
				perform1stConsensusCalling = ((confFile.readline().rstrip()).split("\t"))[1]
				if perform1stConsensusCalling == "yes" or perform1stConsensusCalling == "Yes":
					now = datetime.datetime.now()
					logFile.write("First consensus calling started at "+now.strftime("%H:%M")+"\n\n")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "\nFirst round of consensus calling on sample "+projectName+"\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("mkdir 6_createConsensus")
					os.system("cp 5_refineAssembly/finalScaffold.fasta ./6_createConsensus")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/extractSeqByRange.py ./6_createConsensus")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/revComp.py ./6_createConsensus/")
					os.system("cp "+installationDirectory+"src/scripts/assembly/utils/completeGenome2.py ./6_createConsensus/")
					os.system("cp "+installationDirectory+"src/scripts/utils/biomodule.py ./6_createConsensus/")
					#os.system("cp "+installationDirectory+"src/scripts/utils/biomodule.py ./6_createConsensus/")


					os.chdir("6_createConsensus")

					#Attempts five and three prime ends reconstruction
					os.system(installationDirectory+"src/conda/bin/python completeGenome2.py "+installationDirectory+"  finalScaffold.fasta 0")
					if os.path.isfile("newGenome2.fasta") == True:
						for seq_record in SeqIO.parse("newGenome2.fasta","fasta"):
							newGenome2seq = str(seq_record.seq)
							break
						
						if "N" in newGenome2seq:
							os.system("head -40000 ../1_cleanReads/qualityFiltered_1.fq >subsample_1.fastq")
							os.system("head -40000 ../1_cleanReads/qualityFiltered_2.fq >subsample_2.fastq")
							bowtiePE("newGenome2.fasta","subsample_1.fastq","subsample_2.fastq",self.threadsEntry.get())
							os.system(installationDirectory+"src/conda/bin/picard CollectInsertSizeMetrics I=test_sorted.bam  O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5")
							os.system("head -8 insert_size_metrics.txt | tail -2 | cut -f 6 | tail -1 >insert.size")
							isize = open("insert.size")
							insertSize = isize.readline().rstrip()
							isize.close()
							print(insertSize)
							gfFile = open("gapfillerlib.txt","w")
							gfFile.write("lib1 bwa ../1_cleanReads/qualityFiltered_1.fq ../1_cleanReads/qualityFiltered_2.fq "+ insertSize+" 0.25 FR")
							gfFile.close()
							os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/GapFiller -q "+installationDirectory+" -l gapfillerlib.txt -s newGenome2.fasta -T "+self.threadsEntry.get())
							os.system("cp ./standard_output/standard_output.gapfilled.final.fa newGenome2.fasta")
							for seq_record in SeqIO.parse("newGenome2.fasta","fasta"):
								newGenome2seq = str(seq_record.seq)
								break
						
						os.system("mv finalScaffold.fasta finalScaffold_noEnds.fasta")
						fs = open("finalScaffold.fasta","w")
						fs.write(">finalScaffold\n"+newGenome2seq+"\n")
						fs.close()


					for seq_record in SeqIO.parse("finalScaffold.fasta","fasta"):
						assemblyLength = len(str(seq_record.seq))


					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing first portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python extractSeqByRange.py finalScaffold.fasta finalScaffold 1 15001 f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					bowtiePE("finalScaffold_1_15001_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()

					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")

					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()


					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")# >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_1_15001_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_1_15001_f.txt")

					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_1_15001_f.txt -o output.vcf dedupped.bam")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_1_15001_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_1_15001_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_1_15001_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_firstPortion.vcf")
					os.system("mv output_filtered.vcf  output_firstPortion_filtered.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")


					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing second portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python extractSeqByRange.py finalScaffold.fasta finalScaffold 15001 "+str(assemblyLength -10000 )+" f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					bowtiePE("finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt")
					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt -o output.vcf dedupped.bam")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					#os.system(installationDirectory+"resources/bgzip -c output.vcf_filtered.vcf > output.vcf_filtered.vcf.gz 2>null")
					#os.system(installationDirectory+"resources/tabix output.vcf_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_secondPortion.vcf")
					os.system("mv output_filtered.vcf  output_filtered_secondPortion.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")



					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  Analyzing third portion....\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/python extractSeqByRange.py finalScaffold.fasta finalScaffold "+str(assemblyLength - 10000 )+" 2000000 f")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Aligning reads to the assembly\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					bowtiePE("finalScaffold_"+str(assemblyLength -10000 )+"_2000000_f.txt","../1_cleanReads/qualityFiltered_1.fq","../1_cleanReads/qualityFiltered_2.fq",self.threadsEntry.get())
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Extracting mapped reads\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/samtools view -bF 4 test_sorted.bam >mapped.bam 2>null")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Adding gorup names\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard AddOrReplaceReadGroups I=mapped.bam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Deduplicating\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Calling polymorphisms\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system(installationDirectory+"src/conda/bin/picard CreateSequenceDictionary R=finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/samtools faidx finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt")
					os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.threadsEntry.get()+" -q 30 -Q 30  -f finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt -o output.vcf dedupped.bam")
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/utils/getMajorAllele.py output.vcf output_filtered.vcf >null 2>&1")
					os.system(installationDirectory+"src/conda/bin/perl "+installationDirectory+"src/scripts/assembly/utils/vcf-sort output_filtered.vcf >temp.vcf ; mv temp.vcf output_filtered.vcf")
					os.system(installationDirectory+"src/conda/bin/bgzip -c output_filtered.vcf > output_filtered.vcf.gz 2>null")
					os.system(installationDirectory+"src/conda/bin/tabix output_filtered.vcf.gz >null 2>&1")
					#os.system("java -jar  "+installationDirectory+"resources/GenomeAnalysisTK.jar -T  HaplotypeCaller -R finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt -I dedupped.bam  -o output.vcf -A StrandAlleleCountsBySample >null 2>&1")
					#os.system(installationDirectory+"resources/filterVCF.py output.vcf >null 2>&1")
					#os.system(installationDirectory+"resources/bgzip -c output.vcf_filtered.vcf > output.vcf_filtered.vcf.gz 2>null")
					#os.system(installationDirectory+"resources/tabix output.vcf_filtered.vcf.gz >null 2>&1")
					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "*  *  Creating consensus\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					os.system("cat finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt | "+installationDirectory+"src/conda/bin/bcftools consensus output_filtered.vcf.gz > finalScaffold_"+str(assemblyLength - 10000 )+"_2000000_f.txt_con.fasta 2>null")
					os.system("mv output.vcf output_thirdPortion.vcf")
					os.system("mv output_filtered.vcf  output_filtered_thirdPortion.vcf")
					os.system("rm -f test*")
					os.system("rm -f *.dict")



					finalSequence = ""

					for seq_record in SeqIO.parse("finalScaffold_1_15001_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					for seq_record in SeqIO.parse("finalScaffold_15001_"+str(assemblyLength -10000 )+"_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					for seq_record in SeqIO.parse("finalScaffold_"+str(assemblyLength -10000 )+"_2000000_f.txt_con.fasta","fasta"):
						finalSequence+=str(seq_record.seq)

					outfile = open(projectName+"_genome.fasta","w")
					outfile.write(">finalScaffold\n"+finalSequence+"\n")




					self.logArea.configure(state='normal')
					self.logArea.insert(tk.END, "Second round of consensus calling completed on sample "+projectName+"!!\n")
					self.logArea.see(tk.END)
					self.logArea.configure(state='disabled')
					self.logArea.update()
					now = datetime.datetime.now()
					logFile.write("First consensus calline ended at "+now.strftime("%H:%M")+"\n\n")
					os.chdir("../")
				else:
					if os.path.isfile("./6_createConsensus/"+projectName+"_genome.fasta")==False:
						print("You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......")
						logFile.write("You chose not to run the consensus call step but "+projectName+"_genome.fasta  file is not there. Now exiting......")
						exit()

				self.logArea.configure(state='normal')
				self.logArea.insert(tk.END, "\n\nAll processes completed on sample "+projectName+"!!\n")
				self.logArea.see(tk.END)
				self.logArea.configure(state='disabled')
				self.logArea.update()
				os.chdir("../")



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

		w=670
		h=400
		ws = root.winfo_screenwidth() # width of the screen
		hs = root.winfo_screenheight() # height of the screen
		# calculate x and y coordinates for the Tk root window
		x = (ws/2) - (w/2)
		y = (hs/2) - (h/2)

		top.geometry('%dx%d+%d+%d' % (w, h, x, y))
		top.title("De novo tool")
		top.configure(highlightcolor="black",background="white")

		self.inputFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.inputFileLabel.configure(text="Samples list file")
		self.inputFileLabel.place(x=20,y=20,width=110,height=20)

		self.inputFileEntry = tk.Entry(top)
		self.inputFileEntry.place(x=20,y=40,width=520,height=30)
		self.inputFileEntry.insert(0,"Please select file....")

		self.inputFileButton = tk.Button(top,command=openInputFile,background="#204949",foreground="white")
		self.inputFileButton.place(x=550,y=40,width=100,height=30)
		self.inputFileButton.configure(text="Open file")

		self.runButton = tk.Button(top,command=performAssembly,background="#204949",foreground="white")
		self.runButton.place(x=450,y=350,width=200,height=30)
		self.runButton.configure(text="Run")

		#self.exitButton = tk.Button(top,command=exitProgram)
		#self.exitButton.place(x=550,y=550,width=100,height=30)
		#self.exitButton.configure(text="Exit")

		self.verboseChkValue = tk.BooleanVar()
		self.verboseChkValue.set(False)
		self.verboseCheckButton = tk.Checkbutton(top,variable=self.verboseChkValue,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.verboseCheckButton.place(x=320,y=92,height=20,width=230)
		self.verboseCheckButton.configure(text="Do no delete intermediate files")

		self.threadsEntry = tk.Entry(top)
		self.threadsEntry.place(x=140,y=85,width=30,height=30)
		self.threadsEntry.insert(0,"8")

		self.threadsLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.threadsLabel.configure(text="Number of threads:")
		self.threadsLabel.place(x=20,y=90,width=120,height=20)

		self.memoryLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.memoryLabel.configure(text="Memory(Gb):")
		self.memoryLabel.place(x=190,y=90,width=90,height=20)

		self.memoryEntry = tk.Entry(top)
		self.memoryEntry.place(x=280,y=85,width=30,height=30)
		self.memoryEntry.insert(0,"16")


		self.logFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
		self.logFileLabel.place(x=20,y=120,height=20,width=80)
		self.logFileLabel.configure(text="Log window")

		self.logFrame = tk.Frame(top)
		self.logFrame.place(x=20, y=140, height=240, width=420)
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(borderwidth="2")
		self.logFrame.configure(relief='groove')
		self.logFrame.configure(width=125)
		self.logArea = tk.Text(top,state='disabled')
		self.logArea.place(x=25,y=145,height=230, width=410)
		self.logArea.configure(background="white",borderwidth=5)
		self.logArea.configure(selectbackground="#c4c4c4")


		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Assembly.jpg")
		photo1 =ImageTk.PhotoImage(image)
		logoLabel = tk.Label(top, compound=tk.TOP,height=200,width=200,image=photo1,borderwidth=0,highlightthickness=0)
		logoLabel.place(x=450,y=140)
		logoLabel.image = photo1



if __name__ == '__main__':
	vp_start_gui()
