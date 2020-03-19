# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'readsFiltering.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 

import sys
from os import listdir
from os.path import isfile, join
import os
import time
import numpy as nm
import matplotlib.pyplot as plt

class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Read filtering")
		Form.resize(864, 796)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 20, 261, 20))
		self.label.setObjectName("label")
		self.inputFolderEntry = QtWidgets.QLineEdit(Form)
		self.inputFolderEntry.setGeometry(QtCore.QRect(10, 40, 381, 32))
		self.inputFolderEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.inputFolderEntry.setObjectName("inputFolderEntry")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(10, 90, 261, 20))
		self.label_2.setObjectName("label_2")
		self.outputFolderEntry = QtWidgets.QLineEdit(Form)
		self.outputFolderEntry.setGeometry(QtCore.QRect(10, 110, 381, 32))
		self.outputFolderEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.outputFolderEntry.setObjectName("outputFolderEntry")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(10, 160, 241, 20))
		self.label_3.setObjectName("label_3")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(10, 230, 271, 20))
		self.label_4.setObjectName("label_4")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(10, 300, 241, 20))
		self.label_5.setObjectName("label_5")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(10, 370, 241, 20))
		self.label_6.setObjectName("label_6")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(10, 440, 241, 20))
		self.label_7.setObjectName("label_7")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(20, 510, 241, 20))
		self.label_8.setObjectName("label_8")
		self.recodingFileEntry = QtWidgets.QLineEdit(Form)
		self.recodingFileEntry.setGeometry(QtCore.QRect(10, 180, 381, 32))
		self.recodingFileEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.recodingFileEntry.setObjectName("recodingFileEntry")
		self.bowtieIndexEntry = QtWidgets.QLineEdit(Form)
		self.bowtieIndexEntry.setGeometry(QtCore.QRect(10, 250, 381, 32))
		self.bowtieIndexEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.bowtieIndexEntry.setObjectName("bowtieIndexEntry")
		self.adapter1Entry = QtWidgets.QLineEdit(Form)
		self.adapter1Entry.setGeometry(QtCore.QRect(10, 320, 381, 32))
		self.adapter1Entry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.adapter1Entry.setObjectName("adapter1Entry")
		self.adapter2Entry = QtWidgets.QLineEdit(Form)
		self.adapter2Entry.setGeometry(QtCore.QRect(10, 390, 381, 32))
		self.adapter2Entry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.adapter2Entry.setObjectName("adapter2Entry")
		self.numThreadsEntry = QtWidgets.QLineEdit(Form)
		self.numThreadsEntry.setGeometry(QtCore.QRect(10, 460, 151, 32))
		self.numThreadsEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.numThreadsEntry.setObjectName("numThreadsEntry")
		self.logTextArea = QtWidgets.QTextEdit(Form)
		self.logTextArea.setGeometry(QtCore.QRect(20, 540, 601, 241))
		self.logTextArea.setObjectName("logTextArea")
		self.frame = QtWidgets.QFrame(Form)
		self.frame.setGeometry(QtCore.QRect(560, 40, 291, 231))
		self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
		self.frame.setObjectName("frame")
		self.humanReadsRemovalCheckbox = QtWidgets.QCheckBox(self.frame)
		self.humanReadsRemovalCheckbox.setGeometry(QtCore.QRect(30, 20, 201, 25))
		self.humanReadsRemovalCheckbox.setObjectName("humanReadsRemovalCheckbox")
		self.adaptersTrimmingCheckbox = QtWidgets.QCheckBox(self.frame)
		self.adaptersTrimmingCheckbox.setGeometry(QtCore.QRect(30, 60, 191, 25))
		self.adaptersTrimmingCheckbox.setObjectName("adaptersTrimmingCheckbox")
		self.sampleNameRecodingCheckbox = QtWidgets.QCheckBox(self.frame)
		self.sampleNameRecodingCheckbox.setGeometry(QtCore.QRect(30, 100, 251, 25))
		self.sampleNameRecodingCheckbox.setObjectName("sampleNameRecodingCheckbox")
		self.readsDeduplicationCheckbox = QtWidgets.QCheckBox(self.frame)
		self.readsDeduplicationCheckbox.setGeometry(QtCore.QRect(30, 140, 231, 25))
		self.readsDeduplicationCheckbox.setObjectName("readsDeduplicationCheckbox")
		self.merlinReferenceAlignmentCheckbox = QtWidgets.QCheckBox(self.frame)
		self.merlinReferenceAlignmentCheckbox.setGeometry(QtCore.QRect(30, 180, 251, 25))
		self.merlinReferenceAlignmentCheckbox.setObjectName("merlinReferenceAlignmentCheckbox")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(570, 20, 161, 21))
		self.label_9.setObjectName("label_9")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(650, 540, 201, 211))
		self.label_10.setText("")
		self.label_10.setPixmap(QtGui.QPixmap(installationDirectory+"src/GUI/IconsFinal/filtering.jpg"))
		self.label_10.setObjectName("label_10")
		self.inputFolderButton = QtWidgets.QPushButton(Form)
		self.inputFolderButton.setGeometry(QtCore.QRect(400, 40, 111, 32))
		self.inputFolderButton.setObjectName("inputFolderButton")
		self.outputFolderButton = QtWidgets.QPushButton(Form)
		self.outputFolderButton.setGeometry(QtCore.QRect(400, 110, 112, 32))
		self.outputFolderButton.setObjectName("outputFolderButton")
		self.recodingButton = QtWidgets.QPushButton(Form)
		self.recodingButton.setGeometry(QtCore.QRect(400, 180, 112, 32))
		self.recodingButton.setObjectName("recodingButton")
		self.bowtieIndexButton = QtWidgets.QPushButton(Form)
		self.bowtieIndexButton.setGeometry(QtCore.QRect(400, 250, 112, 32))
		self.bowtieIndexButton.setObjectName("bowtieIndexButton")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(650, 750, 201, 32))
		self.runButton.setObjectName("runButton")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		#Addition to GUI generated by Qt designer
		font = QtGui.QFont()
		font.setPointSize(11)
		Form.setFont(font)
		self.outputFolderButton.clicked.connect(self.selectOutputFolder)
		self.inputFolderButton.clicked.connect(self.selectInputFolder)
		self.recodingButton.clicked.connect(self.selectRecodingFile)
		self.bowtieIndexButton.clicked.connect(self.selectBowtieIndex)
		self.runButton.clicked.connect(lambda:self.runTool(installationDirectory))






	def selectOutputFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select output folder","./")
		self.outputFolderEntry.setText(folderName)

	def selectInputFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select input folder","./")
		self.inputFolderEntry.setText(folderName)

		
	def selectRecodingFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select recoding file","./")
		self.recodingFileEntry.setText(filename)


	def selectBowtieIndex(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select human bowtie2 index","./")
		self.bowtieIndexEntry.setText(filename)

	def plotCoveragePlot(self,covFile,expName):
		position = []
		coverage = []
		windowSize = 0
		avPos = 0
		covFile = open(covFile)
		coverageToPlot = []
		positionToPlot = []
		while True:
			line = covFile.readline().rstrip()
			if not line:
				break
			fields = line.split("\t")
			position.append(int(fields[1]))
			coverage.append(int(fields[2]))
			windowSize +=1
			if windowSize == 200:
				avPos += 200
				positionToPlot.append(avPos)
				avCoverage = nm.mean(coverage)
				coverageToPlot.append(int(avCoverage))
				position = []
				coverage = []
				windowSize = 0

		plt.figure(figsize=(20,10))
		plt.plot(positionToPlot,coverageToPlot,linewidth=0.3)
		plt.xticks(nm.arange(0,235000,10000))
		plt.xlabel("Position (bp)")
		plt.ylabel("Coverage")
		plt.text(1000,950,expName)
		plt.savefig(expName+"_covPlot.png")

	def runTool(self,installationDirectory):
		inputFolder = self.inputFolderEntry.text()
		outputFolder = self.outputFolderEntry.text()
		bowtie2Ref = self.bowtieIndexEntry.text()
		#Collect files to quality check
		onlyfiles = [f for f in listdir(inputFolder) if isfile(join(inputFolder, f))]

		fileNumber = {}
		for item in onlyfiles:
			if not item[:-8] in fileNumber:
				fileNumber[item[:-8]] = 1
			else:
				fileNumber[item[:-8]] += 1

		paired2filter = []
		single2filter = []
		for files in fileNumber:
			if fileNumber[files]==2:
				paired2filter.append(files)
			if fileNumber[files]==1:
				single2filter.append(files)

		#Check recoding file
		codes = {}
		if self.sampleNameRecodingCheckbox.isChecked() == True:
			rcodefile = open(self.recodingFileEntry.text())


			while True:
				line = rcodefile.readline().rstrip()
				if not line:
					break
				fields = line.split("\t")
				oldName = (fields[0].split("/"))[-1]
				if not oldName in codes:
					codes[oldName]=fields[1]
		else:
			for item in paired2filter:
				if not item+"_1.fastq" in codes:
					codes[item+"_1.fastq"] = item+"_1.fastq"
				if not item+"_2.fastq" in codes:
					codes[item+"_2.fastq"] = item+"_2.fastq"



		#Start filtering

		statisticsToReport = []
		statisticsToReport.append("Original reads")
		datasetStatistics = {}
		numTasks = 1
		if self.humanReadsRemovalCheckbox.isChecked() == True:
			numTasks += 4
			statisticsToReport.append("Number Reads after host reads removal")
			statisticsToReport.append("Percentage Reads after host reads removal")
		if self.adaptersTrimmingCheckbox.isChecked() == True:
			numTasks += 2
			statisticsToReport.append("Number Reads after adaptor trimming")
			statisticsToReport.append("Percentage Reads after adaptor trimming")
		if self.readsDeduplicationCheckbox.isChecked() == True:
			numTasks += 2
			statisticsToReport.append("Number Reads after deduplication")
			statisticsToReport.append("Percentage Reads after deduplication")
		if self.merlinReferenceAlignmentCheckbox.isChecked() == True:
			numTasks += 5
			statisticsToReport.append("Merlin coverage for original reads")
			statisticsToReport.append("Number Original reads mapping to Merlin")
			statisticsToReport.append("Percentage Original reads mapping to Merlin")
			statisticsToReport.append("Number of bases covered by original reads")
			statisticsToReport.append("Percentage of bases covered by original reads")


			if self.readsDeduplicationCheckbox.isChecked() == True:
				numTasks += 5
				statisticsToReport.append("Merlin coverage for deduplicated reads")
				statisticsToReport.append("Number of deduplicated reads mapping to Merlin")
				statisticsToReport.append("Percentage of deduplicated reads mapping to Merlin")
				statisticsToReport.append("Number of bases covered by deduplicated reads")
				statisticsToReport.append("Percentage of bases covered by deduplicated reads")




		numTasks = numTasks*len(paired2filter)
		progressBarIncrement = 100000/numTasks
		step = 0

		for dataset in paired2filter:
			if not dataset[0] ==".":
				if not dataset in datasetStatistics:
					datasetStatistics[dataset] = []

				
				self.logTextArea.append("Starting filtering for dataset "+dataset+"")
				self.logTextArea.repaint()
				


				os.system("cp "+inputFolder+"/"+dataset+"_1.fastq tempReads_140875_1.fastq")
				os.system("cp "+inputFolder+"/"+dataset+"_2.fastq tempReads_140875_2.fastq")
				#Check the format of the input fastq file header
				fqfile = open("tempReads_140875_1.fastq")
				header = fqfile.readline().rstrip()
				time.sleep(1)
				if " 1" in header:
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/utils/changeHeaderFormat.py tempReads_140875_1.fastq tempReads_140875_2.fastq")

				
				self.logTextArea.append( "Calculating the number of reads for dataset "+dataset+"")
				self.logTextArea.repaint()
				
				os.system("wc -l "+inputFolder+"/"+dataset+"_1.fastq  >numReads_140875")
				step+=1

				numreads = open("numReads_140875")

				numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
				originalReadsNumber = numberOfReads
				datasetStatistics[dataset].append(numberOfReads)
				numreads.close()
				os.system("rm -f numReads_140875")
				
				self.logTextArea.append( "Done!\nNumber of reads: "+str(numberOfReads)+"")
				self.logTextArea.repaint()
				


				suffixCode = ""
				if self.humanReadsRemovalCheckbox.isChecked() == True:
					suffixCode+="_nh"
					
					self.logTextArea.append( "Removing human reads from dataset "+dataset+"")
					self.logTextArea.repaint()
					
					self.logTextArea.append( "Mapping reads "+dataset+" to the host reference genome....")
					self.logTextArea.repaint()
					os.system(installationDirectory+"src/conda/bin/bowtie2 --local -x "+bowtie2Ref+" -1 tempReads_140875_1.fastq -2 tempReads_140875_2.fastq -p "+self.numThreadsEntry.text()+" -S hostAlignment_140875.sam")

					time.sleep(2)

					
					self.logTextArea.append( "Converting alignment format for reads "+dataset+"")
					self.logTextArea.repaint()
					
					
					os.system(installationDirectory+"src/conda/bin/samtools view -bS -h hostAlignment_140875.sam >hostAlignment_140875.bam")
					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					step+=1

					time.sleep(2)
					
					self.logTextArea.append( "Extracting unmapped reads for dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/bam2fastq --no-aligned --force --strict -o unmapped_140875#.fq hostAlignment_140875.bam")
					step+=1

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					
					os.system("mv unmapped_140875_1.fq tempReads_140875_1.fastq")
					os.system("mv unmapped_140875_2.fq tempReads_140875_2.fastq")
					os.system("rm -f hostAlignment_140875.bam ")

					
					self.logTextArea.append( "Calculating the number of host free reads for dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					os.system("wc -l tempReads_140875_1.fastq >numReads_140875")
					step+=1

					numreads = open("numReads_140875")
					numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
					datasetStatistics[dataset].append(str(numberOfReads))
					datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
					numreads.close()
					os.system("rm -f numReads_140875")
					
					self.logTextArea.append( "Done!...... Number of reads: "+str(numberOfReads)+"")
					self.logTextArea.repaint()
					
					
					



					#step += 1
					#self.progressbar['value']=int( (step/numStep)*100)

					time.sleep(2)

				if self.adaptersTrimmingCheckbox.isChecked() == True:
					suffixCode += "_tr"

					
					self.logTextArea.append( "Trimming reads for dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					print("Performing trim_galore")
					if len(self.adapter1Entry.text())>5 and len(self.adapter2Entry.text())>5:
						adpt1 = self.adapter1Entry.text()
						adpt2 = self.adapter2Entry.text()
						os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt -paired -a "+adpt1+" -a2 "+adpt2+" tempReads_140875_1.fastq tempReads_140875_2.fastq ")

					else:
						os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt -paired tempReads_140875_1.fastq tempReads_140875_2.fastq ")

					step+=1

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					

					os.system("mv tempReads_140875_1_val_1.fq tempReads_140875_1.fastq")
					os.system("mv tempReads_140875_2_val_2.fq tempReads_140875_2.fastq")
					os.system("rm -f tempReads_140875_1.fastq_trimming_report.txt tempReads_140875_2.fastq_trimming_report.txt")

					
					self.logTextArea.append( "Calculating the number reads after trimming for dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					os.system("wc -l tempReads_140875_1.fastq >numReads_140875")
					step+=1

					numreads = open("numReads_140875")
					numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
					datasetStatistics[dataset].append(str(numberOfReads))
					datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
					numreads.close()
					os.system("rm -f numReads_140875")
					
					self.logTextArea.append( "Done!...... Number of reads: "+str(numberOfReads)+"")
					self.logTextArea.repaint()
					
					
					




				if self.readsDeduplicationCheckbox.isChecked() == True:

					
					self.logTextArea.append( "Performing deduplication on dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					os.system("echo tempReads_140875_1.fastq >dedupTemplate.txt")
					os.system("echo tempReads_140875_2.fastq >> dedupTemplate.txt")
					os.system(installationDirectory+"src/conda/bin/fastuniq -i dedupTemplate.txt -t q -o "+codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_dd_1.fastq -p "+codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_dd_2.fastq")
					os.system("rm -f dedupTemplate.txt")
					dedupFileName1 = codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_dd_1.fastq"
					dedupFileName2 = codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_dd_2.fastq"
					
					self.logTextArea.append( "Done!")
					
					
					
					step+=1

					
					self.logTextArea.append( "Calculating the number reads after deduplicating for dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					os.system("wc -l "+dedupFileName1+" >numReads_140875")
					step+=1

					numreads = open("numReads_140875")
					numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
					datasetStatistics[dataset].append(str(numberOfReads))
					datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
					numreads.close()
					os.system("rm -f numReads_140875")
					
					self.logTextArea.append( "Done!...... Number of reads: "+str(numberOfReads)+"")
					
					
					


				if self.merlinReferenceAlignmentCheckbox.isChecked() == True:
					
					self.logTextArea.append( "Performing alignment for original dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/bowtie2 -1 tempReads_140875_1.fastq -2 tempReads_140875_2.fastq -x "+installationDirectory+"data/merlinReference/hcmv -p "+ self.numThreadsEntry.text()+" -S alignment_140875.sam")
					step+=1

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					


					
					self.logTextArea.append( "Convert sam to bam....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools view -bS -h  alignment_140875.sam >alignment_140875.bam")
					step+=1

					self.logTextArea.append( "Done!")
					
					
					


					
					self.logTextArea.append( "Sorting bam....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools sort -o alignment_140875_sorted.bam alignment_140875.bam")
					step+=1
					
					
					
					self.logTextArea.append( "Done!")
					
					
					


					
					self.logTextArea.append( "Calculating coverage for original dataset "+dataset+"....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage_140875.txt")
					os.system(installationDirectory+"src/conda/bin/samtools depth -a -d 10000000 alignment_140875_sorted.bam >coverage_140875.txt")
					self.plotCoveragePlot("coverage_140875.txt",codes[dataset+"_1.fastq"].replace("_1.fastq","")+"_nh_tr")
					os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  | wc -l >breadth_140875.txt")
					step+=1
					
					
					
					self.logTextArea.append( "Done!")
					
					
					


					avCovFile = open("avCoverage_140875.txt")
					avCov = float(avCovFile.readline().rstrip())
					avCovFile.close()
					print("Average coverage",avCov)

					datasetStatistics[dataset].append(avCov)

					os.system(installationDirectory+"src/conda/bin/samtools view -F 4 alignment_140875.sam | wc -l > readsMapping_140875")
					step+=1
					
					
					readsMappingFile = open("readsMapping_140875")
					readsMapping = float(readsMappingFile.readline().rstrip())
					readsMappingFile.close()

					datasetStatistics[dataset].append(str(readsMapping))
					datasetStatistics[dataset].append(str(readsMapping*100/originalReadsNumber)[:4])

					bfile = open("breadth_140875.txt")
					breadthValue = float(bfile.readline().rstrip())
					bfile.close()

					datasetStatistics[dataset].append(str(breadthValue))
					datasetStatistics[dataset].append(str( breadthValue / 235646.0 )[:4])



					if self.readsDeduplicationCheckbox.isChecked() == True:
						
						self.logTextArea.append( "Performing alignment for deduplicated dataset "+dataset+"....")
						self.logTextArea.repaint()
						
						
						os.system(installationDirectory+"src/conda/bin/bowtie2 -1 "+dedupFileName1+" -2 "+ dedupFileName2+ "  -x "+installationDirectory+"data/merlinReference/hcmv -p "+self.numThreadsEntry.text()+" -S alignment_140875.sam")
						step+=1
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Converting sam to bam....")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools view -bS -h  alignment_140875.sam >alignment_140875.bam")
						step+=1
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Sorting bam file....")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools sort -o alignment_140875_sorted.bam alignment_140875.bam")
						step+=1
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Calculating coverage for deduplicated dataset "+dataset+"....")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage_140875.txt")
						os.system(installationDirectory+"src/conda/bin/samtools depth -a -d 10000000 alignment_140875_sorted.bam   >coverage_140875.txt")
						self.plotCoveragePlot("coverage_140875.txt",codes[dataset+"_1.fastq"].replace("_1.fastq","")+"_nh_tr_dd")
						os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  | wc -l >breadth_140875.txt")
						step+=1
						
						
						avCovFile = open("avCoverage_140875.txt")
						avCov = float(avCovFile.readline().rstrip())
						avCovFile.close()
						
						self.logTextArea.append( "Done!")
						self.logTextArea.repaint()
						
						
						

						print("Average coverage for deduplicated ",avCov)

						datasetStatistics[dataset].append(avCov)

						os.system(installationDirectory+"src/conda/bin/samtools view -F 4 alignment_140875.sam | wc -l > readsMapping_140875")
						step+=1
						
						
						readsMappingFile = open("readsMapping_140875")
						readsMapping = float(readsMappingFile.readline().rstrip())
						readsMappingFile.close()

						datasetStatistics[dataset].append(str(readsMapping))
						datasetStatistics[dataset].append(str(readsMapping*100/originalReadsNumber)[:4])

						bfile = open("breadth_140875.txt")
						breadthValue = float(bfile.readline().rstrip())
						bfile.close()

						datasetStatistics[dataset].append(str(breadthValue))
						datasetStatistics[dataset].append(str( breadthValue / 235646.0 )[:4])





				print("finished")
				
				self.logTextArea.append( "\n\nFiltering complete on all datasets!\n")
				self.logTextArea.repaint()
				

				os.system("mv tempReads_140875_1.fastq "+outputFolder+"/"+codes[dataset+"_1.fastq"].replace("_1.fastq","")+suffixCode+"_1.fastq")
				os.system("mv tempReads_140875_2.fastq "+outputFolder+"/"+codes[dataset+"_2.fastq"].replace("_2.fastq","")+suffixCode+"_2.fastq")
				if self.readsDeduplicationCheckbox.isChecked() == True:
					os.system("mv "+dedupFileName1+" "+outputFolder+"/")
					os.system("mv "+dedupFileName2+" "+outputFolder+"/")
				os.system("mv "+codes[dataset+"_1.fastq"].replace("_1.fastq","")+"* "+outputFolder+"/")

				for a in range(len(statisticsToReport)):
					print(statisticsToReport[a],datasetStatistics[dataset][a])

				os.system("rm -f *.sam *.bam readsMapping avCoverage.txt breadth.txt")



				for item in datasetStatistics:
					print(item)
					for a in datasetStatistics[item]:
						print(str(a)+"\t")
				print("\n")


		outfile = open("summaryTable.txt","w")
		outfile.write("Statistics\t")


		for item in datasetStatistics:
			outfile.write(codes[item+"_1.fastq"].replace("_1.fastq","")+"\t")
		outfile.write("\n")
		for a in range(len(statisticsToReport)):
			outfile.write(statisticsToReport[a]+"\t")
			for item in datasetStatistics:
				outfile.write(str(datasetStatistics[item][a])+"\t")
			outfile.write("\n")


		outfile.close()
		os.system("mv summaryTable.txt "+outputFolder+"/")
		os.system("rm -f coverage.txt *140875*")











	











	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Form"))
		self.label.setText(_translate("Form", "Input folder"))
		self.inputFolderEntry.setText(_translate("Form", "No folder selected"))
		self.label_2.setText(_translate("Form", "Output folder"))
		self.outputFolderEntry.setText(_translate("Form", "No folder selected"))
		self.label_3.setText(_translate("Form", "Recoding file"))
		self.label_4.setText(_translate("Form", "Human reference bowtie2 index"))
		self.label_5.setText(_translate("Form", "Adapter 1"))
		self.label_6.setText(_translate("Form", "Adapter 2"))
		self.label_7.setText(_translate("Form", "Number of threads"))
		self.label_8.setText(_translate("Form", "Log window"))
		self.recodingFileEntry.setText(_translate("Form", "No file selected"))
		self.bowtieIndexEntry.setText(_translate("Form", "No index selected"))
		self.adapter1Entry.setText(_translate("Form", "Auto detect"))
		self.adapter2Entry.setText(_translate("Form", "Auto detect"))
		self.numThreadsEntry.setText(_translate("Form", "8"))
		self.humanReadsRemovalCheckbox.setText(_translate("Form", "Human reads removal"))
		self.adaptersTrimmingCheckbox.setText(_translate("Form", "Adapters trimming"))
		self.sampleNameRecodingCheckbox.setText(_translate("Form", "Sample file name recoding"))
		self.readsDeduplicationCheckbox.setText(_translate("Form", "Reads deduplication"))
		self.merlinReferenceAlignmentCheckbox.setText(_translate("Form", "Merlin reference alignment"))
		self.label_9.setText(_translate("Form", "Tasks to perform"))
		self.inputFolderButton.setText(_translate("Form", "Open folder"))
		self.outputFolderButton.setText(_translate("Form", "Open folder"))
		self.recodingButton.setText(_translate("Form", "Open file"))
		self.bowtieIndexButton.setText(_translate("Form", "Open file"))
		self.runButton.setText(_translate("Form", "Run"))


if __name__ == "__main__":
	import sys

	installationDirectory = sys.argv[1]
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())

