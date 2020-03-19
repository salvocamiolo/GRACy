# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dbsubmission.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog 
from PyQt5.QtWidgets import QInputDialog
import os,sys,time

class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(848, 623)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(460, 20, 281, 20))
		self.label.setObjectName("label")
		self.sampleToSubmitEntry = QtWidgets.QLineEdit(Form)
		self.sampleToSubmitEntry.setGeometry(QtCore.QRect(460, 40, 241, 32))
		self.sampleToSubmitEntry.setObjectName("sampleToSubmitEntry")
		self.sampleToSubmitButton = QtWidgets.QPushButton(Form)
		self.sampleToSubmitButton.setGeometry(QtCore.QRect(710, 40, 112, 32))
		self.sampleToSubmitButton.setObjectName("sampleToSubmitButton")
		self.inputFolderButton = QtWidgets.QPushButton(Form)
		self.inputFolderButton.setGeometry(QtCore.QRect(270, 40, 112, 32))
		self.inputFolderButton.setObjectName("inputFolderButton")
		self.label_2 = QtWidgets.QLabel(Form)
		self.label_2.setGeometry(QtCore.QRect(20, 20, 281, 20))
		self.label_2.setObjectName("label_2")
		self.inputFolderEntry = QtWidgets.QLineEdit(Form)
		self.inputFolderEntry.setGeometry(QtCore.QRect(20, 40, 241, 32))
		self.inputFolderEntry.setObjectName("inputFolderEntry")
		self.projectInfoButton = QtWidgets.QPushButton(Form)
		self.projectInfoButton.setGeometry(QtCore.QRect(270, 120, 112, 32))
		self.projectInfoButton.setObjectName("projectInfoButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(20, 100, 281, 20))
		self.label_3.setObjectName("label_3")
		self.projectInfoEntry = QtWidgets.QLineEdit(Form)
		self.projectInfoEntry.setGeometry(QtCore.QRect(20, 120, 241, 32))
		self.projectInfoEntry.setObjectName("projectInfoEntry")
		self.sampleInfoButton = QtWidgets.QPushButton(Form)
		self.sampleInfoButton.setGeometry(QtCore.QRect(710, 120, 112, 32))
		self.sampleInfoButton.setObjectName("sampleInfoButton")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(460, 100, 281, 20))
		self.label_4.setObjectName("label_4")
		self.sampleInfoEntry = QtWidgets.QLineEdit(Form)
		self.sampleInfoEntry.setGeometry(QtCore.QRect(460, 120, 241, 32))
		self.sampleInfoEntry.setObjectName("sampleInfoEntry")
		self.readsInfoButton = QtWidgets.QPushButton(Form)
		self.readsInfoButton.setGeometry(QtCore.QRect(270, 200, 112, 32))
		self.readsInfoButton.setObjectName("readsInfoButton")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(20, 180, 281, 20))
		self.label_5.setObjectName("label_5")
		self.readsInfoEntry = QtWidgets.QLineEdit(Form)
		self.readsInfoEntry.setGeometry(QtCore.QRect(20, 200, 241, 32))
		self.readsInfoEntry.setObjectName("readsInfoEntry")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(460, 180, 141, 20))
		self.label_6.setObjectName("label_6")
		self.usernameEntry = QtWidgets.QLineEdit(Form)
		self.usernameEntry.setGeometry(QtCore.QRect(460, 200, 171, 32))
		self.usernameEntry.setObjectName("usernameEntry")
		self.passwordEntry = QtWidgets.QLineEdit(Form)
		self.passwordEntry.setGeometry(QtCore.QRect(650, 200, 171, 32))
		self.passwordEntry.setEchoMode(QtWidgets.QLineEdit.Password)
		self.passwordEntry.setObjectName("passwordEntry")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(650, 180, 141, 20))
		self.label_7.setObjectName("label_7")
		self.newProjectCombo = QtWidgets.QComboBox(Form)
		self.newProjectCombo.setGeometry(QtCore.QRect(20, 280, 90, 28))
		self.newProjectCombo.setObjectName("newProjectCombo")
		self.newProjectCombo.addItem("")
		self.newProjectCombo.addItem("")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(20, 260, 161, 20))
		self.label_8.setObjectName("label_8")
		self.newSampleCombo = QtWidgets.QComboBox(Form)
		self.newSampleCombo.setGeometry(QtCore.QRect(230, 280, 90, 28))
		self.newSampleCombo.setObjectName("newSampleCombo")
		self.newSampleCombo.addItem("")
		self.newSampleCombo.addItem("")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(230, 260, 181, 20))
		self.label_9.setObjectName("label_9")
		self.testSubmissionCombo = QtWidgets.QComboBox(Form)
		self.testSubmissionCombo.setGeometry(QtCore.QRect(460, 280, 90, 28))
		self.testSubmissionCombo.setObjectName("testSubmissionCombo")
		self.testSubmissionCombo.addItem("")
		self.testSubmissionCombo.addItem("")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(460, 260, 161, 20))
		self.label_10.setObjectName("label_10")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(20, 340, 76, 20))
		self.label_11.setObjectName("label_11")
		self.logArea = QtWidgets.QTextEdit(Form)
		self.logArea.setGeometry(QtCore.QRect(20, 360, 581, 241))
		self.logArea.setObjectName("logArea")
		self.label_12 = QtWidgets.QLabel(Form)
		self.label_12.setGeometry(QtCore.QRect(620, 350, 211, 221))
		self.label_12.setText("")
		self.label_12.setPixmap(QtGui.QPixmap("GRACy_easyinstall/src/GUI/IconsFinal/DB.jpg"))
		self.label_12.setObjectName("label_12")
		self.pushButton_6 = QtWidgets.QPushButton(Form)
		self.pushButton_6.setGeometry(QtCore.QRect(620, 570, 201, 32))
		self.pushButton_6.setObjectName("pushButton_6")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)


		#Addition to GUI generated by Qt designer
		self.inputFolderButton.clicked.connect(self.selectInputFolder)
		self.sampleToSubmitButton.clicked.connect(self.selectSampleToSubmit)
		self.projectInfoButton.clicked.connect(self.selectProjectTable)
		self.sampleInfoButton.clicked.connect(self.selectSampleTable)
		self.readsInfoButton.clicked.connect(self.selectReadsInfo)
		self.sampleToSubmitButton.clicked.connect(self.submitSamples)

	def selectInputFolder(self):
		foldername = QFileDialog.getExistingDirectory(None,"Select input folder","./")
		self.inputFolderEntry.setText(foldername)
	
	def selectSampleToSubmit(self):
		filename,__ = QFileDialog.getOpenFileName(None,"Select File","./")
		self.sampleToSubmitEntry.setText(filename)

	def selectProjectTable(self):
		filename,__ = QFileDialog.getOpenFileName(None,"Select project info file","./")
		self.projectInfoEntry.setText(filename)

	def selectSampleTable(self):
		filename,__ = QFileDialog.getOpenFileName(None,"Select sample info file","./")
		self.sampleInfoEntry.setText(filename)

	def selectReadsInfo(self):
		filename,__ = QFileDialog.getOpenFileName(None,"Select reads info file","./")
		self.readsInfoEntry.setText(filename)

	def createNewProject(self,projectFile):
		if os.path.isfile(projectFile)==True:
			projectInfoFile = open(projectFile)
			projectInfoFile.readline()
			line = projectInfoFile.readline().rstrip()
			projectInfoFields = line.split("\t")
			projectAlias = projectInfoFields[0]
			projectName = projectInfoFields[1]
			projectTitle = projectInfoFields[2]
			projectDescription = projectInfoFields[3]
			project_xml = open("project.xml","w")
			project_xml.write("<PROJECT_SET>\n")
			project_xml.write("<PROJECT alias=\""+projectAlias+"\">\n")
			project_xml.write("<NAME>"+projectName+"</NAME>\n")
			project_xml.write("<TITLE>"+projectTitle+"</TITLE>\n")
			project_xml.write("<DESCRIPTION>"+projectDescription+"</DESCRIPTION>\n")
			project_xml.write("<SUBMISSION_PROJECT>\n<SEQUENCING_PROJECT></SEQUENCING_PROJECT>\n</SUBMISSION_PROJECT>\n<PROJECT_LINKS>\n<PROJECT_LINK>\n<XREF_LINK>\n<DB></DB>\n<ID></ID>\n</XREF_LINK>\n</PROJECT_LINK>\n</PROJECT_LINKS>\n</PROJECT>\n</PROJECT_SET>")
			#IMPORTANT! To change username and password in GRACy
			project_xml.close()
			if self.testSubmissionchkValue.text()==True:
				os.system("curl -u "+self.usernameEntry.text()+":"+self.passwordEntry.text()+" -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt")
			else:
				MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
				if MsgBox == 'yes':
					os.system("curl -u "+self.usernameEntry.text()+":"+self.passwordEntry.text()+" -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt")
				else:
					exit()
			receiptFile = open("projectReceipt")

			while True:
				line = receiptFile.readline().rstrip()
				if not line:
					break
				fields = line.split("\"")
				for item in fields:
					if "PRJE" in item:
						receiptFile.close()
						return item
						break

	def createNewSample(self,sampleName):
		if not sampleName in sampleAccession:
			sampleAccession[sampleName] = ""
		infile = open(self.sampleInfoEntry.text())

		infile.readline() #Read header
		
		sampleFound = 0
		while True:
			line = infile.readline().rstrip()
			#print line
			#sys.stdin.read(1)
			if not line:
				break
			fields = line.split("\t")
			alias = fields[0]
			center = fields[1]
			title = fields[2]
			taxonID = fields[3]
			scientificName = fields[4]
			geographicLocation = fields[5]
			hostCommonName = fields[6]
			hostHealthState = fields[7]
			isolationSource = fields[8]
			sex = fields[9]
			hostScientificName = fields[10]
			collectorName = fields[11]
			collectingInstitution = fields[12]
			isolate = fields[13]
			
			if alias == sampleName:
				sampleFound = 1
				xmlfile = open(alias+"_sample.xml","w")
				xmlfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<SAMPLE_SET>\n<SAMPLE alias=\""+alias+"\" center_name=\""+center+"\">\n")
				xmlfile.write("<TITLE>"+title+"</TITLE>\n<SAMPLE_NAME>\n<TAXON_ID>"+taxonID+"</TAXON_ID>\n<SCIENTIFIC_NAME>"+scientificName+"</SCIENTIFIC_NAME>\n")
				xmlfile.write(" <COMMON_NAME></COMMON_NAME>\n</SAMPLE_NAME>\n<SAMPLE_ATTRIBUTES>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write(" <TAG>geographic location (country and/or sea)</TAG>\n<VALUE>"+geographicLocation+"</VALUE>\n")
				xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>host common name</TAG>\n<VALUE>"+hostCommonName+"</VALUE>\n")
				xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>host health state</TAG>\n<VALUE>"+hostHealthState+"</VALUE>\n")
				xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>isolation source host associated</TAG>\n<VALUE>"+isolationSource+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write("<TAG>host sex</TAG>\n<VALUE>"+sex+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write("<TAG>host scientific name</TAG>\n<VALUE>"+hostScientificName+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write("<TAG>collector name</TAG>\n<VALUE>"+collectorName+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write("<TAG>collecting institution</TAG>\n<VALUE>"+collectingInstitution+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
				xmlfile.write("<TAG>isolate</TAG>\n<VALUE>"+isolate+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n</SAMPLE_ATTRIBUTES>\n</SAMPLE>\n</SAMPLE_SET>\n")
				xmlfile.close()
				if self.testSubmissionchkValue.text()==True:
					os.system("curl -u "+self.usernameEntry.text()+":"+self.passwordEntry.text()+" -F \"SUBMISSION=@submission.xml\" -F \"SAMPLE=@"+alias+"_sample.xml\" \"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/\" >sampleReceipt.txt")
				else:
					MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
					if MsgBox == 'yes':
						os.system("curl -u "+self.usernameEntry.text()+":"+self.passwordEntry.text()+" -F \"SUBMISSION=@submission.xml\" -F \"SAMPLE=@"+alias+"_sample.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >sampleReceipt.txt")

				receiptFile = open("sampleReceipt.txt")
				while True:
					line = receiptFile.readline().rstrip()
					if not line:
						break
					fields = line.split("\"")
					for item in fields:
						if "ERS" in item or "SRS" in item or "DRS" in item:
							receiptFile.close()
							return item
							break

		if sampleFound == 0:
			print("Sample",sampleName,"was not found in the provided sample file.")
			print ("Now exit....")
			exit()

		infile.close()


	def submitSamples(self):
		if self.sampleToSubmitEntry.text()=="Please select a file...." or self.inputFolderEntry.text()=="Please select a folder...." or self.readsInfoEntry.text()=="Please select a file...." or self.usernameEntry.text()=="myENA_Username" or self.passwordEntry.text()=="":
			
			#self.logArea.append("The following fields must be filled:\n")
			#self.logArea.append("Samples to submit\nFastq folder\nFastq info folder\nENA username\nPassword\n\n")
			exit()
			
			
		else:
			#*******************************************************************
			#********************* Main algorithm reads submission  ************
			#*******************************************************************
			isTest = self.testSubmissionCombo.currentText()
			inputFolder = self.inputFolderEntry.text()
			filesToSubmit = self.sampleToSubmitEntry.text()
			createProject = self.newProjectCombo.currentText()
			createSample = self.newSampleCombo.currentText()
			projectInfo = self.projectInfoEntry.text()
			sampleInfo = self.sampleInfoEntry.text()
			readsInfo = self.readsInfoEntry.text()
			os.system("cp "+installationDirectory+"src/scripts/dbsubmission/utils/submission.xml .")

			

			
			self.logArea.append("isTest: "+str(isTest)+"\n"+"inputFolder: "+str(inputFolder)+"\nfilesToSubmit: "+str(filesToSubmit)+"\nCreatePoject: "+str(createProject)+"\nCreateSample: "+str(createSample)+"\nProjectInfo: "+str(projectInfo)+"\nSampleInfo: "+str(sampleInfo)+"\nReadsInfo: "+str(readsInfo)+"\n\n")
			self.logArea.repaint()
			
			

			if createProject=="Yes":
				
				self.logArea.append("Creating a new Project....\n")
				self.logArea.repaint()
				
				
				
				projectAccessionNumber = createNewProject(projectInfo)
				
				self.logArea.append("A project was created on ENA with the accession number: "+projectAccessionNumber+"\n\n")
				self.logArea.repaint()
				
				
				

			toSubmitFile = open(filesToSubmit)
			while True:
				sampleName = toSubmitFile.readline().rstrip()
				if not sampleName:
					break

				if not sampleName in experimentAccession:
					experimentAccession[sampleName] = ""

				if not sampleName in runAccession:
					runAccession[sampleName] = ""

				sampleList.append(sampleName)


				#Create samples if needed
				if createSample == "Yes":
					
					self.logArea.append("Creating a new sample for dataset "+sampleName+"....\n")
					self.logArea.repaint()
					
					
					
					sampleAccessionNumber = createNewSample(sampleName)
					
					self.logArea.append("Sample created for dataset "+sampleName+" with accession number: "+sampleAccessionNumber+"\n")
					self.logArea.repaint()
					
					
					

					time.sleep(3)
				
				#Grab reads information from readsInfo file
				infile = open(readsInfo)
				infile.readline()
				foundRecord = 0
				while True:
					line = infile.readline().rstrip()
					if not line:
						break
					fields = line.split("\t")
					projectAccessionField = fields[0]
					sampleAccessionField = fields[1]

					sampleAlias = fields[2]
					instrument = fields[3]
					insertSize = fields[4]
					librarySource = fields[5]
					librarySelection = fields[6]
					libraryStrategy = fields[7]
					fq1 = fields[8]
					fq2 = fields[9]

					if createProject == "Yes":
						projectAccessionField = projectAccessionNumber
					if createSample =="Yes":
						sampleAccessionField = sampleAccessionNumber

					sampleAccession[sampleName] = sampleAccessionField

					if sampleAlias == sampleName:
						foundRecord = 1
						os.system("ln -s "+inputFolder+"/"+fq1)
						os.system("ln -s "+inputFolder+"/"+fq2)

						manifest = open(sampleName+"_manifestFile.txt","w",buffering=0)
						manifest.write("INSTRUMENT\t"+instrument+"\nINSERT_SIZE\t"+insertSize+"\nLIBRARY_SOURCE\t"+librarySource+"\nLIBRARY_SELECTION\t"+librarySelection+"\n")
						manifest.write("STUDY\t"+projectAccessionField+"\nSAMPLE\t"+sampleAccessionField+"\nNAME\t"+sampleName+"\n")
						manifest.write("LIBRARY_STRATEGY\t"+libraryStrategy+"\nFASTQ\t"+fq1+"\nFASTQ\t"+fq2+"\n")
						manifest.close()
						
						
						self.logArea.append("Submitting reads for sample "+sampleName+"....\n")
						self.logArea.repaint()
						
						
						

						if self.testSubmissionCombo.currentText()=="Yes":
							os.system(installationDirectory+"src/conda/bin/java -Xmx2048m -jar "+installationDirectory+"src/scripts/dbsubmission/utils/webin-cli-2.2.0.jar  -context reads -userName "+self.usernameEntry.text()+" -password "+self.passwordEntry.text()+"  -manifest "+sampleName+"_manifestFile.txt -test -submit >fastqReceipt")
						else:
							MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
							if MsgBox == 'yes':
								os.system(installationDirectory+"src/conda/bin/java -Xmx2048m -jar "+installationDirectory+"src/scripts/dbsubmission/utils/webin-cli-2.2.0.jar -context reads -userName "+self.usernameEntry.text()+" -password "+self.passwordEntry.text()+"  -manifest "+sampleName+"_manifestFile.txt -submit >fastqReceipt")
							else:
								exit()
						#os.system("rm -f "+fq1+" "+fq2)
						#print "Waiting ENA output...."
						#time.sleep(10)
						os.system("grep . fastqReceipt > tempF ; mv tempF fastqReceipt")
						receiptFile = open("fastqReceipt")
						while True:
							line = receiptFile.readline().rstrip()
							if not line:
								break
							fields = line.split(" ")
							print(fields)
							for item in fields:
								if "ERR" in item or "SRR" in item or "DRR" in item: 
									if not "ERROR" in item:
										runAccession[sampleName] = item

								if "ERX" in item or "SRX" in item or "DRX" in item: 
									experimentAccession[sampleName]=item

								
						receiptFile.close()
						
						self.logArea.append("Fastq files submitted!\n")
						self.logArea.append("Run accession number: "+runAccession[sampleName]+"\n")
						self.logArea.append("Experiment accession number: "+experimentAccession[sampleName]+"\n")
						self.logArea.repaint()
						
						
						
						
						print(experimentAccession[sampleName],"Submitted")
						sys.stdin.read(1)



				print("Finished")
				print(sampleName,sampleAccession[sampleName],runAccession[sampleName],experimentAccession[sampleName])





			outfile = open("accessionNumbersTable.txt","w")
			outfile.write("SampleAlias\tSampleAccession\tRunAccession\tExperimentAccession\n")

		for item in sampleList:
			#print(item,sampleAccession[item],runAccession[item],experimentAccession[item])
			outfile.write(item+"\t"+sampleAccession[item]+"\t"+runAccession[item]+"\t"+experimentAccession[item]+"\n")

		outfile.close()


	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Sequencing data submission tool"))
		self.label.setText(_translate("Form", "Samples to submit"))
		self.sampleToSubmitButton.setText(_translate("Form", "Open file"))
		self.inputFolderButton.setText(_translate("Form", "Open file"))
		self.label_2.setText(_translate("Form", "Inpupt folder"))
		self.projectInfoButton.setText(_translate("Form", "Open file"))
		self.label_3.setText(_translate("Form", "Project info table"))
		self.sampleInfoButton.setText(_translate("Form", "Open file"))
		self.label_4.setText(_translate("Form", "Samples info table"))
		self.readsInfoButton.setText(_translate("Form", "Open file"))
		self.label_5.setText(_translate("Form", "Read datasets info table"))
		self.label_6.setText(_translate("Form", "ENA username"))
		self.label_7.setText(_translate("Form", "ENA password"))
		self.newProjectCombo.setItemText(0, _translate("Form", "Yes"))
		self.newProjectCombo.setItemText(1, _translate("Form", "No"))
		self.label_8.setText(_translate("Form", "Create new project"))
		self.newSampleCombo.setItemText(0, _translate("Form", "Yes"))
		self.newSampleCombo.setItemText(1, _translate("Form", "No"))
		self.label_9.setText(_translate("Form", "Create new sample(s)"))
		self.testSubmissionCombo.setItemText(0, _translate("Form", "Yes"))
		self.testSubmissionCombo.setItemText(1, _translate("Form", "No"))
		self.label_10.setText(_translate("Form", "Test submission"))
		self.label_11.setText(_translate("Form", "Log area"))
		self.pushButton_6.setText(_translate("Form", "Submit"))


if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)
	Form = QtWidgets.QWidget()
	installationDirectory = sys.argv[1]
	ui = Ui_Form()
	ui.setupUi(Form,installationDirectory)
	Form.show()
	sys.exit(app.exec_())

