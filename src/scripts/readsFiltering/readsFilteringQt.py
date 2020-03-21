# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'readsFiltering.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will process lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QMessageBox

import sys
from os import listdir
from os.path import isfile, join
import os
import time
import numpy as nm
import matplotlib.pyplot as plt



class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Reads filtering tool")
		Form.resize(1160, 718)
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(440, 10, 261, 20))
		self.label_2.setObjectName("label_2")
		self.outputFolderEntry = QtWidgets.QLineEdit(Form)
		self.outputFolderEntry.setGeometry(QtCore.QRect(440, 30, 251, 32))
		self.outputFolderEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.outputFolderEntry.setObjectName("outputFolderEntry")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(440, 80, 241, 20))
		self.label_3.setObjectName("label_3")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(440, 150, 271, 20))
		self.label_4.setObjectName("label_4")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(440, 220, 241, 20))
		self.label_5.setObjectName("label_5")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(440, 290, 241, 20))
		self.label_6.setObjectName("label_6")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(440, 360, 241, 20))
		self.label_7.setObjectName("label_7")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(20, 430, 241, 20))
		self.label_8.setObjectName("label_8")
		self.recodingFileEntry = QtWidgets.QLineEdit(Form)
		self.recodingFileEntry.setGeometry(QtCore.QRect(440, 100, 251, 32))
		self.recodingFileEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.recodingFileEntry.setObjectName("recodingFileEntry")
		self.bowtieIndexEntry = QtWidgets.QLineEdit(Form)
		self.bowtieIndexEntry.setGeometry(QtCore.QRect(440, 170, 251, 32))
		self.bowtieIndexEntry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.bowtieIndexEntry.setObjectName("bowtieIndexEntry")
		self.adapter1Entry = QtWidgets.QLineEdit(Form)
		self.adapter1Entry.setGeometry(QtCore.QRect(440, 240, 251, 32))
		self.adapter1Entry.setText("")
		self.adapter1Entry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.adapter1Entry.setObjectName("adapter1Entry")
		self.adapter2Entry = QtWidgets.QLineEdit(Form)
		self.adapter2Entry.setGeometry(QtCore.QRect(440, 310, 251, 32))
		self.adapter2Entry.setText("")
		self.adapter2Entry.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
		self.adapter2Entry.setObjectName("adapter2Entry")
		self.logTextArea = QtWidgets.QTextEdit(Form)
		self.logTextArea.setGeometry(QtCore.QRect(20, 460, 901, 241))
		self.logTextArea.setObjectName("logTextArea")
		self.frame = QtWidgets.QFrame(Form)
		self.frame.setGeometry(QtCore.QRect(850, 30, 291, 231))
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
		self.label_9.setGeometry(QtCore.QRect(860, 10, 161, 21))
		self.label_9.setObjectName("label_9")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(940, 460, 201, 211))
		self.label_10.setText("")
		self.label_10.setPixmap(QtGui.QPixmap(installationDirectory+"src/GUI/IconsFinal/filtering.jpg"))
		self.label_10.setObjectName("label_10")
		self.outputFolderButton = QtWidgets.QPushButton(Form)
		self.outputFolderButton.setGeometry(QtCore.QRect(700, 30, 112, 32))
		self.outputFolderButton.setObjectName("outputFolderButton")
		self.recodingButton = QtWidgets.QPushButton(Form)
		self.recodingButton.setGeometry(QtCore.QRect(700, 100, 112, 32))
		self.recodingButton.setObjectName("recodingButton")
		self.bowtieIndexButton = QtWidgets.QPushButton(Form)
		self.bowtieIndexButton.setGeometry(QtCore.QRect(700, 170, 112, 32))
		self.bowtieIndexButton.setObjectName("bowtieIndexButton")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(940, 670, 201, 32))
		self.runButton.setObjectName("runButton")
		self.numThreadsCombo = QtWidgets.QComboBox(Form)
		self.numThreadsCombo.setGeometry(QtCore.QRect(440, 380, 90, 28))
		self.numThreadsCombo.setObjectName("numThreadsCombo")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.numThreadsCombo.addItem("")
		self.selectedFilesArea = QtWidgets.QTextEdit(Form)
		self.selectedFilesArea.setGeometry(QtCore.QRect(20, 30, 271, 381))
		self.selectedFilesArea.setObjectName("selectedFilesArea")
		self.selectFilesButton = QtWidgets.QPushButton(Form)
		self.selectFilesButton.setGeometry(QtCore.QRect(300, 30, 112, 32))
		self.selectFilesButton.setObjectName("selectFilesButton")
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(20, 10, 101, 20))
		self.label.setObjectName("label")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		#Addition to GUI generated by Qt designer

		self.outputFolderButton.clicked.connect(self.selectOutputFolder)
		self.selectFilesButton.clicked.connect(self.selectFiles)
		self.recodingButton.clicked.connect(self.selectRecodingFile)
		self.bowtieIndexButton.clicked.connect(self.selectBowtieIndex)
		self.runButton.clicked.connect(lambda:self.runTool(installationDirectory))
		self.recodingFileEntry.setReadOnly(True)
		self.outputFolderEntry.setReadOnly(True)
		self.bowtieIndexEntry.setReadOnly(True)







	onlyfiles = []
	def selectOutputFolder(self):
		folderName = QFileDialog.getExistingDirectory(None, "Select output folder","./")
		if len(folderName)>2:
			self.outputFolderEntry.setText(folderName)

	def selectFiles(self):
		self.selectedFilesArea.clear()
		self.onlyfiles = []
		filenames,__ = QFileDialog.getOpenFileNames(None, "Select paired end fastq files","./")
		filenames = sorted(filenames)
		if len(filenames)>0:
			for a in filenames:
				if "_1.fastq" in a:
					self.selectedFilesArea.append((a.replace("_1.fastq","").split("/")[-1]))
				if "R1_001.fastq" in a:
					self.selectedFilesArea.append((a.replace("_R1_001.fastq","").split("/"))[-1])
				if "_1.fq" in a:
					self.selectedFilesArea.append((a.replace("_1.fq","").split("/"))[-1])
				if "R1_001.fq" in a:
					self.selectedFilesArea.append((a.replace("_R1_001.fq","").split("/"))[-1])

		
		
		for a in range(0,len(filenames)-1,+2):
			if ("_1.fastq" in filenames[a] or "_2.fastq" in filenames[a] or "_R1_001.fastq" in filenames[a] or "_R2_001.fastq" in filenames[a] or "_1.fq" in filenames[a] or "_2.fq" in filenames[a] or "_R1_001.fq" in filenames[a] or "_R2_001.fq" in filenames[a])  and ( "_1.fastq" in filenames[a+1] or "_2.fastq" in filenames[a+1] or "_R1_001.fastq" in filenames[a+1] or "_R2_001.fastq" in filenames[a+1] or "_1.fq" in filenames[a+1] or "_2.fq" in filenames[a+1] or "_R1_001.fq" in filenames[a+1] or "_R2_001.fq" in filenames[a+1]):
				self.onlyfiles.append((filenames[a],filenames[a+1]))
			else:
				msg = QMessageBox()
				msg.setIcon(QMessageBox.Warning)
				msg.setText("Some of the selected files are not in the expected format")
				msg.setWindowTitle("Warning")
				msg.setDetailedText("Accepted format are _1.fastq   _2.fastq\n_1.fq   _2.fastq\n_R1_001.fastq   _R2_001.fastq\n_R1_001.fq   _R2_001.fq\n ")
				msg.setStandardButtons(QMessageBox.Ok)
				msg.exec_()


	def refreshTextArea(self,selected):
		self.selectedFilesArea.clear()
		for item in self.onlyfiles:
			if "_1.fastq" in item[0]:

				if (item[0].replace("_1.fastq","").split("/"))[-1] == selected:
					self.selectedFilesArea.append((item[0].replace("_1.fastq","").split("/")[-1])+"  <--- ")
				else:
					self.selectedFilesArea.append((item[0].replace("_1.fastq","").split("/")[-1]))
			if "R1_001.fastq" in item[0]:
				if (item[0].replace("_R1_001.fastq","").split("/"))[-1] == selected:
					self.selectedFilesArea.append((item[0].replace("_R1_001.fastq","").split("/"))[-1]+"  <---")
				else:
					self.selectedFilesArea.append((item[0].replace("_R1_001.fastq","").split("/"))[-1])
			if "_1.fq" in item[0]:
				if (item[0].replace("_1.fq","").split("/"))[-1] == selected:
					self.selectedFilesArea.append((item[0].replace("_1.fq","").split("/"))[-1]+"  <---")
				else:
					self.selectedFilesArea.append((item[0].replace("_1.fq","").split("/"))[-1])
			if "R1_001.fq" in item[0]:
				if (item[0].replace("_R1_001.fq","").split("/"))[-1] == selected:
					self.selectedFilesArea.append((item[0].replace("_R1_001.fq","").split("/"))[-1]+"  <---")
				else:
					self.selectedFilesArea.append((item[0].replace("_R1_001.fq","").split("/"))[-1])
			
		
		


	def getPrefix(self,a):
		if "_1.fastq" in a:
			return ((a.replace("_1.fastq","").split("/")[-1]))
		if "R1_001.fastq" in a:
			return ((a.replace("_R1_001.fastq","").split("/"))[-1])
		if "_1.fq" in a:
			return ((a.replace("_1.fq","").split("/"))[-1])
		if "R1_001.fq" in a:
			return ((a.replace("_R1_001.fq","").split("/"))[-1])



		
	def selectRecodingFile(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select recoding file","./")
		if len(filename)>2:
			self.recodingFileEntry.setText(filename)


	def selectBowtieIndex(self):
		filename, __ = QFileDialog.getOpenFileName(None,"Select human bowtie2 index","./","*index.1.bt2")

		if len(filename)>2:
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
		plt.title(expName)
		plt.savefig(expName+"_covPlot.png")

	def runTool(self,installationDirectory):

	
		if str(self.selectedFilesArea.toPlainText) == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("No fastq files has been selected.")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select at least one paired end reads dataset.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		outputFolder = self.outputFolderEntry.text()
		if outputFolder == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("An output folder should be selected.")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should open the folder where all the produced output files will be saved")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		recodFile = self.recodingFileEntry.text()
		if recodFile == "" and self.sampleNameRecodingCheckbox.isChecked()==True:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid recoding table should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("It seems like you chose to recode the your sequencing data filenames but you did not provide any recoding table file.\n You should select a file reporting two tab separated columns (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		dnaBases = "agctAGCT"
		adapter1Chars = self.adapter1Entry.text()
		adapter2Chars = self.adapter2Entry.text()
		for base in dnaBases:
			adapter1Chars = adapter1Chars.replace(base,"")
			adapter2Chars = adapter2Chars.replace(base,"")
		if len(adapter1Chars)>0 or len(adapter2Chars) >0:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("One of the adapter contains non-canonical DNA bases")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("Only bases A, G, C and T are allowed in the adapter string ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return


		bowtie2Ref = self.bowtieIndexEntry.text()
		if bowtie2Ref == "" and self.humanReadsRemovalCheckbox.isChecked()==True:
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid bowtie2 index for the human reference genome should be specified")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("It seems like you chose to remove the human reads but no bowtie2 index for the human reference genome was specified. Please select the file ending with the suffix \"index.1.bt2\" ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return


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
			for item in self.onlyfiles:
				justFile = item[0].split("/")[-1]
				if not justFile in codes:
					codes[justFile] = justFile
				justFile = item[1].split("/")[-1]
				if not justFile in codes:
					codes[justFile] = justFile


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




		for dataset in self.onlyfiles:
			if not dataset[0] ==".":
				if not dataset in datasetStatistics:
					datasetStatistics[dataset] = []
				
				self.refreshTextArea(self.getPrefix( (dataset[0].split("/"))[-1]))
				for item in self.onlyfiles:
					print(item)

				self.logTextArea.append("Starting filtering for reads in dataset "+self.getPrefix(codes[(dataset[0].split("/"))[-1]]))
				self.logTextArea.repaint()
				


				os.system("cp "+dataset[0] + " tempReads_140875_1.fastq")
				os.system("cp "+dataset[1] + " tempReads_140875_2.fastq")
				#Check the format of the input fastq file header
				fqfile = open("tempReads_140875_1.fastq")
				header = fqfile.readline().rstrip()
				time.sleep(1)
				if " 1" in header:
					os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/utils/changeHeaderFormat.py tempReads_140875_1.fastq tempReads_140875_2.fastq")

				
				self.logTextArea.append( "Calculating the number of reads")
				self.logTextArea.repaint()
				
				os.system("wc -l "+dataset[0]+"   >numReads_140875")
				

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
					
					self.logTextArea.append( "Removing human reads from dataset "+self.getPrefix(codes[(dataset[0].split("/"))[-1]]))
					self.logTextArea.repaint()
					
					self.logTextArea.append( "Mapping reads to the host reference genome....")
					self.logTextArea.repaint()


					os.system(installationDirectory+"src/conda/bin/bowtie2 --local -x "+bowtie2Ref.replace(".1.bt2","")+" -1 tempReads_140875_1.fastq -2 tempReads_140875_2.fastq -p "+self.numThreadsCombo.currentText()+" -S hostAlignment_140875.sam")
					#subprocess.call(installationDirectory+"src/conda/bin/bowtie2 --local -x "+bowtie2Ref.replace(".1.bt2","")+" -1 tempReads_140875_1.fastq -2 tempReads_140875_2.fastq -p "+self.numThreadsCombo.currentText()+" -S hostAlignment_140875.sam",shell=True)

					time.sleep(2)

					
					self.logTextArea.append( "Converting alignment format for reads")
					self.logTextArea.repaint()
					
					
					os.system(installationDirectory+"src/conda/bin/samtools view -bS -h hostAlignment_140875.sam >hostAlignment_140875.bam")
					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					

					time.sleep(2)
					
					self.logTextArea.append( "Extracting unmapped reads")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/bam2fastq --no-aligned --force --strict -o unmapped_140875#.fq hostAlignment_140875.bam")
					

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					
					os.system("mv unmapped_140875_1.fq tempReads_140875_1.fastq")
					os.system("mv unmapped_140875_2.fq tempReads_140875_2.fastq")
					os.system("rm -f hostAlignment_140875.bam ")

					
					self.logTextArea.append( "Calculating the number of host free reads")
					self.logTextArea.repaint()
					
					
					os.system("wc -l tempReads_140875_1.fastq >numReads_140875")
					

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

					
					self.logTextArea.append( "Trimming reads for dataset "+self.getPrefix(codes[(dataset[0].split("/"))[-1]]))
					self.logTextArea.repaint()
					
					
					
					print("Performing trim_galore")
					if len(self.adapter1Entry.text())>5 and len(self.adapter2Entry.text())>5:
						adpt1 = self.adapter1Entry.text()
						adpt2 = self.adapter2Entry.text()
						os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt -paired -a "+adpt1+" -a2 "+adpt2+" tempReads_140875_1.fastq tempReads_140875_2.fastq ")

					else:
						os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt -paired tempReads_140875_1.fastq tempReads_140875_2.fastq ")

					

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					

					os.system("mv tempReads_140875_1_val_1.fq tempReads_140875_1.fastq")
					os.system("mv tempReads_140875_2_val_2.fq tempReads_140875_2.fastq")
					os.system("rm -f tempReads_140875_1.fastq_trimming_report.txt tempReads_140875_2.fastq_trimming_report.txt")

					
					self.logTextArea.append( "Calculating the number reads after trimming")
					self.logTextArea.repaint()
					
					
					os.system("wc -l tempReads_140875_1.fastq >numReads_140875")
					

					numreads = open("numReads_140875")
					numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
					datasetStatistics[dataset].append(str(numberOfReads))
					datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
					numreads.close()
					os.system("rm -f numReads_140875")
					
					self.logTextArea.append( "Done!...... Number of reads: "+str(numberOfReads)+"")
					self.logTextArea.repaint()
					
					
					




				if self.readsDeduplicationCheckbox.isChecked() == True:

					
					self.logTextArea.append( "Performing deduplication on dataset "+self.getPrefix(codes[(dataset[0].split("/"))[-1]]))
					self.logTextArea.repaint()
					
					
					
					os.system("echo tempReads_140875_1.fastq >dedupTemplate.txt")
					os.system("echo tempReads_140875_2.fastq >> dedupTemplate.txt")
					os.system(installationDirectory+"src/conda/bin/fastuniq -i dedupTemplate.txt -t q -o "+ self.getPrefix( codes[ (dataset[0].split("/"))[-1] ] )+suffixCode+"_dd_1.fastq -p "+self.getPrefix( codes[ (dataset[0].split("/"))[-1] ]) +suffixCode+"_dd_2.fastq")
					os.system("rm -f dedupTemplate.txt")
					dedupFileName1 = self.getPrefix( codes[ (dataset[0].split("/"))[-1] ] )+suffixCode+"_dd_1.fastq"
					dedupFileName2 = self.getPrefix( codes[ (dataset[0].split("/"))[-1] ] )+suffixCode+"_dd_2.fastq"
					
					self.logTextArea.append( "Done!")
					
					
					
					

					
					self.logTextArea.append( "Calculating the number reads after deduplicating for dataset "+self.getPrefix(codes[(dataset[0].split("/"))[-1]]))
					self.logTextArea.repaint()
					
					
					
					os.system("wc -l "+dedupFileName1+" >numReads_140875")
					

					numreads = open("numReads_140875")
					numberOfReads = int(((numreads.readline().rstrip()).split(" "))[0])/2
					datasetStatistics[dataset].append(str(numberOfReads))
					datasetStatistics[dataset].append(str(numberOfReads*100/originalReadsNumber)[:4])
					numreads.close()
					os.system("rm -f numReads_140875")
					
					self.logTextArea.append( "Done!...... Number of reads: "+str(numberOfReads)+"")
					
					
					


				if self.merlinReferenceAlignmentCheckbox.isChecked() == True:
					
					self.logTextArea.append( "Performing alignment for original dataset")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/bowtie2 -1 tempReads_140875_1.fastq -2 tempReads_140875_2.fastq -x "+installationDirectory+"data/merlinReference/hcmv -p "+ self.numThreadsCombo.currentText()+" -S alignment_140875.sam")
					

					
					self.logTextArea.append( "Done!")
					self.logTextArea.repaint()
					
					
					


					
					self.logTextArea.append( "Convert sam to bam....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools view -bS -h  alignment_140875.sam >alignment_140875.bam")
					

					self.logTextArea.append( "Done!")
					
					
					


					
					self.logTextArea.append( "Sorting bam....")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools sort -o alignment_140875_sorted.bam alignment_140875.bam")
					
					
					
					
					self.logTextArea.append( "Done!")
					
					
					


					
					self.logTextArea.append( "Calculating coverage for original dataset ")
					self.logTextArea.repaint()
					
					
					
					os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage_140875.txt")
					os.system(installationDirectory+"src/conda/bin/samtools depth -a -d 10000000 alignment_140875_sorted.bam >coverage_140875.txt")
					self.plotCoveragePlot("coverage_140875.txt",self.getPrefix(codes[  (dataset[0].split("/"))[-1] ] ))
					os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  | wc -l >breadth_140875.txt")
					
					
					
					
					self.logTextArea.append( "Done!")
					
					
					


					avCovFile = open("avCoverage_140875.txt")
					avCov = float(avCovFile.readline().rstrip())
					avCovFile.close()
					print("Average coverage",avCov)

					datasetStatistics[dataset].append(avCov)

					os.system(installationDirectory+"src/conda/bin/samtools view -F 4 alignment_140875.sam | wc -l > readsMapping_140875")
					
					
					
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
						
						self.logTextArea.append( "Performing alignment for deduplicated dataset")
						self.logTextArea.repaint()
						
						
						os.system(installationDirectory+"src/conda/bin/bowtie2 -1 "+dedupFileName1+" -2 "+ dedupFileName2+ "  -x "+installationDirectory+"data/merlinReference/hcmv -p "+self.numThreadsCombo.currentText()+" -S alignment_140875.sam")
						
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Converting sam to bam....")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools view -bS -h  alignment_140875.sam >alignment_140875.bam")
						
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Sorting bam file....")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools sort -o alignment_140875_sorted.bam alignment_140875.bam")
						
						
						
						
						self.logTextArea.append( "Done!")
						
						
						

						
						self.logTextArea.append( "Calculating coverage for deduplicated dataset ")
						self.logTextArea.repaint()
						
						
						
						os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  |  awk '{sum+=$3} END { print sum/NR}' >avCoverage_140875.txt")
						os.system(installationDirectory+"src/conda/bin/samtools depth -a -d 10000000 alignment_140875_sorted.bam   >coverage_140875.txt")
						self.plotCoveragePlot("coverage_140875.txt",self.getPrefix(codes[(dataset[0].split("/"))[-1]])+"_nh_tr_dd")
						os.system(installationDirectory+"src/conda/bin/samtools depth -d 10000000 alignment_140875_sorted.bam  | wc -l >breadth_140875.txt")
						
						
						
						avCovFile = open("avCoverage_140875.txt")
						avCov = float(avCovFile.readline().rstrip())
						avCovFile.close()
						
						self.logTextArea.append( "Done!")
						self.logTextArea.repaint()
						
						
						

						print("Average coverage for deduplicated ",avCov)

						datasetStatistics[dataset].append(avCov)

						os.system(installationDirectory+"src/conda/bin/samtools view -F 4 alignment_140875.sam | wc -l > readsMapping_140875")
						
						
						
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


				os.system("mv tempReads_140875_1.fastq "+outputFolder+"/"+self.getPrefix(codes[(dataset[0].split("/"))[-1]])  +suffixCode+"_1.fastq")
				os.system("mv tempReads_140875_2.fastq "+outputFolder+"/"+self.getPrefix(codes[(dataset[0].split("/"))[-1]]) +suffixCode+"_2.fastq")
				if self.readsDeduplicationCheckbox.isChecked() == True:
					os.system("mv "+dedupFileName1+" "+outputFolder+"/")
					os.system("mv "+dedupFileName2+" "+outputFolder+"/")
				os.system("mv "+self.getPrefix(codes[(dataset[0].split("/"))[-1]])+"* "+outputFolder+"/")

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
			outfile.write(self.getPrefix(codes[(dataset[0].split("/"))[-1]])  +"\t")
		outfile.write("\n")
		for a in range(len(statisticsToReport)):
			outfile.write(statisticsToReport[a]+"\t")
			for item in datasetStatistics:
				outfile.write(str(datasetStatistics[item][a])+"\t")
			outfile.write("\n")


		outfile.close()
		os.system("mv summaryTable.txt "+outputFolder+"/")
		os.system("rm -f coverage.txt *140875*")
		print("finished")
				
		self.logTextArea.append( "\n\nFiltering complete on all datasets!\n")
		self.logTextArea.repaint()


		

	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Reads filtering tool"))
		self.label_2.setText(_translate("Form", "Output folder"))
		self.outputFolderEntry.setText(_translate("Form", ""))
		self.label_3.setText(_translate("Form", "Recoding file"))
		self.label_4.setText(_translate("Form", "Human reference bowtie2 index"))
		self.label_5.setText(_translate("Form", "Adapter 1"))
		self.label_6.setText(_translate("Form", "Adapter 2"))
		self.label_7.setText(_translate("Form", "Number of threads"))
		self.label_8.setText(_translate("Form", "Log window"))
		self.recodingFileEntry.setText(_translate("Form", ""))
		self.bowtieIndexEntry.setText(_translate("Form", ""))
		self.humanReadsRemovalCheckbox.setText(_translate("Form", "Human reads removal"))
		self.adaptersTrimmingCheckbox.setText(_translate("Form", "Adapters trimming"))
		self.sampleNameRecodingCheckbox.setText(_translate("Form", "Sample file name recoding"))
		self.readsDeduplicationCheckbox.setText(_translate("Form", "Reads deduplication"))
		self.merlinReferenceAlignmentCheckbox.setText(_translate("Form", "Merlin reference alignment"))
		self.label_9.setText(_translate("Form", "Tasks to perform"))
		self.outputFolderButton.setText(_translate("Form", "Open folder"))
		self.recodingButton.setText(_translate("Form", "Open file"))
		self.bowtieIndexButton.setText(_translate("Form", "Open file"))
		self.runButton.setText(_translate("Form", "Run"))
		self.numThreadsCombo.setItemText(0, _translate("Form", "1"))
		self.numThreadsCombo.setItemText(1, _translate("Form", "2"))
		self.numThreadsCombo.setItemText(2, _translate("Form", "3"))
		self.numThreadsCombo.setItemText(3, _translate("Form", "4"))
		self.numThreadsCombo.setItemText(4, _translate("Form", "5"))
		self.numThreadsCombo.setItemText(5, _translate("Form", "6"))
		self.numThreadsCombo.setItemText(6, _translate("Form", "7"))
		self.numThreadsCombo.setItemText(7, _translate("Form", "8"))
		self.numThreadsCombo.setItemText(8, _translate("Form", "9"))
		self.numThreadsCombo.setItemText(9, _translate("Form", "10"))
		self.numThreadsCombo.setItemText(10, _translate("Form", "11"))
		self.numThreadsCombo.setItemText(11, _translate("Form", "12"))
		self.numThreadsCombo.setItemText(12, _translate("Form", "13"))
		self.numThreadsCombo.setItemText(13, _translate("Form", "14"))
		self.numThreadsCombo.setItemText(14, _translate("Form", "15"))
		self.numThreadsCombo.setItemText(15, _translate("Form", "16"))
		self.numThreadsCombo.setItemText(16, _translate("Form", "17"))
		self.numThreadsCombo.setItemText(17, _translate("Form", "18"))
		self.numThreadsCombo.setItemText(18, _translate("Form", "19"))
		self.numThreadsCombo.setItemText(19, _translate("Form", "20"))
		self.numThreadsCombo.setItemText(20, _translate("Form", "21"))
		self.numThreadsCombo.setItemText(21, _translate("Form", "22"))
		self.numThreadsCombo.setItemText(22, _translate("Form", "23"))
		self.numThreadsCombo.setItemText(23, _translate("Form", "24"))
		self.selectFilesButton.setText(_translate("Form", "Select files"))
		self.label.setText(_translate("Form", "Fastq files"))


if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	installationDirectory = sys.argv[1]
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())

