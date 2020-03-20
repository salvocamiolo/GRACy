# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../Dropbox/Projects/GRACy/GRACy_UI/snpCalling.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QInputDialog
import sys,os

class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(1222, 580)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(10, 10, 211, 20))
		self.label.setObjectName("label")
		self.inputFileButton = QtWidgets.QPushButton(Form)
		self.inputFileButton.setGeometry(QtCore.QRect(460, 30, 112, 32))
		self.inputFileButton.setObjectName("inputFileButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(600, 80, 191, 20))
		self.label_2.setObjectName("label_2")
		self.cdsFileEntry = QtWidgets.QLineEdit(Form)
		self.cdsFileEntry.setGeometry(QtCore.QRect(600, 100, 241, 32))
		self.cdsFileEntry.setObjectName("cdsFileEntry")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(600, 220, 191, 20))
		self.label_3.setObjectName("label_3")
		self.annotationFileEntry = QtWidgets.QLineEdit(Form)
		self.annotationFileEntry.setGeometry(QtCore.QRect(600, 240, 241, 32))
		self.annotationFileEntry.setObjectName("annotationFileEntry")
		self.cdsFileButton = QtWidgets.QPushButton(Form)
		self.cdsFileButton.setGeometry(QtCore.QRect(850, 100, 112, 32))
		self.cdsFileButton.setObjectName("cdsFileButton")
		self.annotationFileButton = QtWidgets.QPushButton(Form)
		self.annotationFileButton.setGeometry(QtCore.QRect(850, 240, 112, 32))
		self.annotationFileButton.setObjectName("annotationFileButton")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(600, 150, 191, 20))
		self.label_4.setObjectName("label_4")
		self.genomeFileEntry = QtWidgets.QLineEdit(Form)
		self.genomeFileEntry.setGeometry(QtCore.QRect(600, 170, 241, 32))
		self.genomeFileEntry.setObjectName("genomeFileEntry")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(600, 10, 191, 20))
		self.label_5.setObjectName("label_5")
		self.genomeFileButton = QtWidgets.QPushButton(Form)
		self.genomeFileButton.setGeometry(QtCore.QRect(850, 170, 112, 32))
		self.genomeFileButton.setObjectName("genomeFileButton")
		self.outputFolderEntry = QtWidgets.QLineEdit(Form)
		self.outputFolderEntry.setGeometry(QtCore.QRect(600, 30, 241, 32))
		self.outputFolderEntry.setObjectName("outputFolderEntry")
		self.outputFolderButton = QtWidgets.QPushButton(Form)
		self.outputFolderButton.setGeometry(QtCore.QRect(850, 30, 112, 32))
		self.outputFolderButton.setObjectName("outputFolderButton")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(1050, 10, 191, 20))
		self.label_6.setObjectName("label_6")
		self.numThreadsCombo = QtWidgets.QComboBox(Form)
		self.numThreadsCombo.setGeometry(QtCore.QRect(1110, 30, 90, 28))
		self.numThreadsCombo.setLayoutDirection(QtCore.Qt.LeftToRight)
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
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(1050, 80, 191, 20))
		self.label_7.setObjectName("label_7")
		self.baseQualityCombo = QtWidgets.QComboBox(Form)
		self.baseQualityCombo.setGeometry(QtCore.QRect(1110, 100, 90, 28))
		self.baseQualityCombo.setObjectName("baseQualityCombo")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.baseQualityCombo.addItem("")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(10, 300, 191, 20))
		self.label_8.setObjectName("label_8")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(1000, 320, 201, 211))
		self.label_9.setText("")
		self.label_9.setPixmap(QtGui.QPixmap(installationDirectory+"src/GUI/IconsFinal/SNPcalling.jpg"))
		self.label_9.setObjectName("label_9")
		self.logArea = QtWidgets.QTextEdit(Form)
		self.logArea.setGeometry(QtCore.QRect(10, 320, 981, 241))
		self.logArea.setObjectName("logArea")
		self.runButton = QtWidgets.QPushButton(Form)
		self.runButton.setGeometry(QtCore.QRect(1000, 530, 201, 32))
		self.runButton.setObjectName("runButton")
		self.deduplicationCombo = QtWidgets.QComboBox(Form)
		self.deduplicationCombo.setGeometry(QtCore.QRect(1110, 170, 90, 28))
		self.deduplicationCombo.setObjectName("deduplicationCombo")
		self.deduplicationCombo.addItem("")
		self.deduplicationCombo.addItem("")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(1020, 150, 191, 20))
		self.label_10.setObjectName("label_10")
		self.selectedFilesArea = QtWidgets.QTextEdit(Form)
		self.selectedFilesArea.setGeometry(QtCore.QRect(10, 30, 441, 251))
		self.selectedFilesArea.setObjectName("selectedFilesArea")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		#Addition to GUI generated by Qt designer

		self.inputFileButton.clicked.connect(self.selectFiles)
		self.cdsFileButton.clicked.connect(self.selectCDSfile)
		self.annotationFileButton.clicked.connect(self.selectAnnotationFile)
		self.genomeFileButton.clicked.connect(self.selectGenomeFile)
		self.outputFolderButton.clicked.connect(self.selectOuputFolder)
		self.runButton.clicked.connect(self.runTool)


	onlyfiles = []
	def selectFiles(self):
		self.selectedFilesArea.clear()
		self.onlyfiles = []
		filenames,__ = QFileDialog.getOpenFileNames(None, "Select paired end fastq files","./")
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

	
	def selectCDSfile(self):
		filename,__ = QFileDialog.getOpenFileName(None,"Select CDS file","./")
		self.cdsFileEntry.setText(filename) 

	def selectAnnotationFile(self):
		filename, __= QFileDialog.getOpenFileName(None,"Select annotation file","./")
		self.annotationFileEntry.setText(filename)
	
	def selectGenomeFile(self):
		filename, __= QFileDialog.getOpenFileName(None,"Select Genome fasta file","./")
		self.genomeFileEntry.setText(filename)

	def selectOuputFolder(self):
		foldername = QFileDialog.getExistingDirectory(None,"Select output folder","./")
		self.outputFolderEntry.setText(foldername)

	def runTool(self):
		outputFolder = self.outputFolderEntry.text()
		if outputFolder == "No folder selected":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("An output folder should be selected.")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should open the folder where all the produced output files will be saved")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		if str(self.selectedFilesArea.toPlainText) == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("No fastq files has been selected.")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select at least one paired end reads dataset.\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		referenceFile = self.genomeFileEntry.text()
		if referenceFile == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A reference genome should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should open the fasta formatted sequence file of the genome you are working on")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		gffFile = self.annotationFileEntry.text()
		if gffFile == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("An annotation file should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should open the gff formatted annotation file of the genome you are working on")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		cdsFile = self.cdsFileEntry.text()
		if cdsFile == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A CDS file should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should open the fasta formatted coding sequence file of the genome you are working on")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return
			
		

				

		for dataset in self.onlyfiles:
			read1 = dataset[0]
			read2 = dataset[1]
			sampleName = self.getPrefix((dataset[0].split("/"))[-1])
			self.refreshTextArea(sampleName)
			

			sampleTempFolder = sampleName+"_4864985345_tempFolder"
			os.system("mkdir -p "+sampleTempFolder)
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/conda/bin/bowtie2-build "+referenceFile+" ./"+sampleTempFolder+"/reference -q")

			print("Analising sample",sampleName)
			
			self.logArea.append("Analising sample"+sampleName+"....")
			self.logArea.repaint()
			
			
			


			
			self.logArea.append("*  Reads quality filtering before alignment")
			self.logArea.repaint()
			
			
			

			
			self.logArea.append("*  *  Running Trimgalore....")
			self.logArea.repaint()
			
			
			



			print("Reads quality filtering before alignment")
			print("Running Trimgalore")
			if ".fastq" in read1:
				prefix1 = ((read1.split("/"))[-1]).replace(".fastq","")
				prefix2 = ((read2.split("/"))[-1]).replace(".fastq","")
			if ".fq" in read1:
				prefix1 = ((read1.split("/"))[-1]).replace(".fq","")
				prefix2 = ((read2.split("/"))[-1]).replace(".fq","")

			os.system(installationDirectory+"src/conda/bin/trim_galore --path_to_cutadapt "+installationDirectory+"src/conda/bin/cutadapt --paired -q 0  "+read1+" "+read2 + " -o "+sampleTempFolder)

			
			self.logArea.append("*  *  Performing deduplication....")
			self.logArea.repaint()
			
			
			


			print("Performing deduplication")
			os.system("echo "+sampleTempFolder+"/"+prefix1+"_val_1.fq > "+sampleTempFolder+"/inputFastUniq")
			os.system("echo "+sampleTempFolder+"/"+prefix2+"_val_2.fq >> "+sampleTempFolder+"/inputFastUniq")
			os.system(installationDirectory+"src/conda/bin/fastuniq -i "+sampleTempFolder+"/inputFastUniq -t q -o "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_1.fastq -p "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_2.fastq")

			print("Performing prinseq quality filtering")

			
			self.logArea.append("*  Running Prinseq....")
			self.logArea.repaint()
			
			
			
			os.system(installationDirectory+"src/conda/bin/prinseq-lite.pl -fastq "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_1.fastq  -fastq2 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_2.fastq -min_qual_mean 25 -trim_qual_right 30 -trim_ns_right 20  -trim_qual_window 5 -trim_qual_step 1 -min_len 80 -out_bad null -out_good "+sampleName+"_trimmed_dedup_pr")
			os.system("mv "+sampleName+"_trimmed_dedup_pr* ./"+sampleTempFolder)
			os.system(installationDirectory+"/src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/trimPolyN.py "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_1.fastq "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_2.fastq")

			
			self.logArea.append("*  Reads alignment on reference")
			self.logArea.repaint()
			
			
			

			
			self.logArea.append("*  *  Performing alignment....")
			self.logArea.repaint()
			
			
			

			print("Aigning reads to reference")
			os.system(installationDirectory+"src/conda/bin/bowtie2 --end-to-end -1 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_1.fastq -2 "+sampleTempFolder+"/"+sampleName+"_trimmed_dedup_pr_2.fastq -x "+sampleTempFolder+"/reference -S "+sampleTempFolder+"/"+sampleName+"_alignment.sam -p "+self.numThreadsCombo.currentText())

			
			self.logArea.append("*  *  Converting sam to bam....")
			self.logArea.repaint()
			
			
			

			print("Converting sam to bam")
			os.system(installationDirectory+"src/conda/bin/samtools view -bS -h "+sampleTempFolder+"/"+sampleName+"_alignment.sam > "+sampleTempFolder+"/"+sampleName+"_alignment.bam")

			
			self.logArea.append("*  *  Sorting bam....")
			self.logArea.repaint()
			
			
			

			print("Sorting bam")
			os.system(installationDirectory+"src/conda/bin/samtools sort -o "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+sampleTempFolder+"/"+sampleName+"_alignment.bam")

			if self.deduplicationCombo.currentText()=="Yes":
				
				self.logArea.append("*  *  Marking duplicates....")
				self.logArea.repaint()
				
				
				

				print("Marking duplicates ")
				os.system(installationDirectory+"src/conda/bin/picard  AddOrReplaceReadGroups I="+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam O="+sampleTempFolder+"/"+sampleName+"_rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")
				os.system(installationDirectory+"src/conda/bin/picard MarkDuplicates I="+sampleTempFolder+"/"+sampleName+"_rg_added_sorted.bam O="+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M="+sampleTempFolder+"/output.metrics")

			else:
				os.system("mv "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam")
			print("Calling snps with lofreq")
			#Analyze snps with lowfreq
			os.system(installationDirectory+"/src/conda/bin/samtools faidx "+referenceFile)
			
			self.logArea.append("*  Calling SNPs with lofreq")
			self.logArea.repaint()
			
			
			

			os.system(installationDirectory+"src/conda/bin/lofreq call-parallel --pp-threads "+self.numThreadsCombo.currentText()+" -q "+str(self.baseQualityCombo.currentText())  +" -Q "+str(self.baseQualityCombo.currentText()) +" -f "+referenceFile+" -o "+sampleTempFolder+"/"+sampleName+"_SNPs.vcf "+sampleTempFolder+"/"+sampleName+"_markedDuplicates.bam")
			print("Lofreq terminated")
			#sys.stdin.read(1)

			#Analyze indes using the GATK pipeline

			


			os.system(installationDirectory+"src/conda/bin/samtools faidx "+referenceFile)




			os.system("mv "+sampleTempFolder+"/"+sampleName+"_alignment_sorted.bam "+outputFolder+"/")
			os.system("rm "+sampleTempFolder+"/"+sampleName+"_alignment* -f")
			#os.system("rm -f inputFastUniq "+prefix1+"_val_1.fq "+prefix2+"_val_2.fq  "+sampleName+"_trimmed_dedup_* "+sampleName+"_markedDuplicates.ba* "+sampleName+"_rg_added_sorted.bam "+sampleName+"*vcf.idx*")
			os.system(" ls "+sampleTempFolder+"/*SNPs.vcf > "+sampleTempFolder+"/vcfFilesToExamine")
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/polyAn.py "+sampleTempFolder+"/vcfFilesToExamine"+" "+self.cdsFileEntry.text()+" "+self.annotationFileEntry.text())
			#os.system("rm vcfFilesToExamine *.dict *.fai *trimming_report.txt")
			os.system("ls "+sampleTempFolder+"/*_snpEffect.txt > "+sampleTempFolder+"/file2Plot")
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/buildTable.py "+sampleTempFolder+"/file2Plot")
			os.system("mv  "+sampleTempFolder+"/*_snpEffect.txt "+sampleTempFolder+"/*_snpFreq.txt "+sampleTempFolder+"/*.vcf "+outputFolder+"/" )
			os.system("mv "+sampleTempFolder+"/* ./"+outputFolder+"/")
			os.system("rm -rf "+sampleTempFolder)

		
		
		#os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/polyAn.py "+self.inputFileEntry.text()+" "+self.cdsFileEntry.text()+" "+self.annotationFileEntry.text())
		#os.system("rm vcfFilesToExamine *.dict  mapped.bam  rg_added_sorted.bam *trimming_report.txt")
		#os.system("ls *_snpEffect.txt > file2Plot")
		#os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/utils/buildTable.py file2Plot")
		#os.system("mv  *_snpEffect.txt *_snpFreq.txt snpTable.txt "+outputFolder+"/" )




	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "SNP calling tool"))
		self.label.setText(_translate("Form", "Input file"))
		self.inputFileButton.setText(_translate("Form", "Select files"))
		self.label_2.setText(_translate("Form", "CDS fasta file"))
		self.label_3.setText(_translate("Form", "Annotation gff file"))
		self.cdsFileButton.setText(_translate("Form", "Open file"))
		self.annotationFileButton.setText(_translate("Form", "Open file"))
		self.label_4.setText(_translate("Form", "Genome fasta file"))
		self.label_5.setText(_translate("Form", "Output folder"))
		self.genomeFileButton.setText(_translate("Form", "Open file"))
		self.outputFolderButton.setText(_translate("Form", "Open file"))
		self.label_6.setText(_translate("Form", "Number of threads"))
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
		self.label_7.setText(_translate("Form", "Base quality cutoff"))
		self.baseQualityCombo.setItemText(0, _translate("Form", "30"))
		self.baseQualityCombo.setItemText(1, _translate("Form", "31"))
		self.baseQualityCombo.setItemText(2, _translate("Form", "32"))
		self.baseQualityCombo.setItemText(3, _translate("Form", "33"))
		self.baseQualityCombo.setItemText(4, _translate("Form", "34"))
		self.baseQualityCombo.setItemText(5, _translate("Form", "35"))
		self.baseQualityCombo.setItemText(6, _translate("Form", "36"))
		self.baseQualityCombo.setItemText(7, _translate("Form", "37"))
		self.baseQualityCombo.setItemText(8, _translate("Form", "38"))
		self.baseQualityCombo.setItemText(9, _translate("Form", "39"))
		self.baseQualityCombo.setItemText(10, _translate("Form", "40"))
		self.baseQualityCombo.setItemText(11, _translate("Form", "29"))
		self.baseQualityCombo.setItemText(12, _translate("Form", "28"))
		self.baseQualityCombo.setItemText(13, _translate("Form", "27"))
		self.baseQualityCombo.setItemText(14, _translate("Form", "26"))
		self.baseQualityCombo.setItemText(15, _translate("Form", "25"))
		self.baseQualityCombo.setItemText(16, _translate("Form", "24"))
		self.baseQualityCombo.setItemText(17, _translate("Form", "23"))
		self.baseQualityCombo.setItemText(18, _translate("Form", "22"))
		self.baseQualityCombo.setItemText(19, _translate("Form", "21"))
		self.baseQualityCombo.setItemText(20, _translate("Form", "20"))
		self.baseQualityCombo.setItemText(21, _translate("Form", "19"))
		self.baseQualityCombo.setItemText(22, _translate("Form", "18"))
		self.baseQualityCombo.setItemText(23, _translate("Form", "17"))
		self.baseQualityCombo.setItemText(24, _translate("Form", "16"))
		self.baseQualityCombo.setItemText(25, _translate("Form", "15"))
		self.baseQualityCombo.setItemText(26, _translate("Form", "14"))
		self.baseQualityCombo.setItemText(27, _translate("Form", "13"))
		self.baseQualityCombo.setItemText(28, _translate("Form", "12"))
		self.baseQualityCombo.setItemText(29, _translate("Form", "11"))
		self.baseQualityCombo.setItemText(30, _translate("Form", "10"))
		self.baseQualityCombo.setItemText(31, _translate("Form", "9"))
		self.baseQualityCombo.setItemText(32, _translate("Form", "8"))
		self.baseQualityCombo.setItemText(33, _translate("Form", "7"))
		self.baseQualityCombo.setItemText(34, _translate("Form", "6"))
		self.baseQualityCombo.setItemText(35, _translate("Form", "5"))
		self.baseQualityCombo.setItemText(36, _translate("Form", "4"))
		self.baseQualityCombo.setItemText(37, _translate("Form", "3"))
		self.baseQualityCombo.setItemText(38, _translate("Form", "2"))
		self.baseQualityCombo.setItemText(39, _translate("Form", "1"))
		self.baseQualityCombo.setItemText(40, _translate("Form", "0"))
		self.label_8.setText(_translate("Form", "Log area"))
		self.runButton.setText(_translate("Form", "Run"))
		self.deduplicationCombo.setItemText(0, _translate("Form", "Yes"))
		self.deduplicationCombo.setItemText(1, _translate("Form", "No"))
		self.label_10.setText(_translate("Form", "Perform deduplication"))


if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	installationDirectory = sys.argv[1]
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())

