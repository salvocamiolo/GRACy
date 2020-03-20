# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../Dropbox/Projects/GRACy/GRACy_UI/dbsubmission.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
	def setupUi(self, Form,installationDirectory):
		Form.setObjectName("Form")
		Form.resize(1163, 619)
		self.label = QtWidgets.QLabel(Form)
		self.label.setGeometry(QtCore.QRect(20, 20, 281, 20))
		self.label.setObjectName("label")
		self.sampleToSubmitButton = QtWidgets.QPushButton(Form)
		self.sampleToSubmitButton.setGeometry(QtCore.QRect(330, 40, 112, 32))
		self.sampleToSubmitButton.setObjectName("sampleToSubmitButton")
		self.projectInfoButton = QtWidgets.QPushButton(Form)
		self.projectInfoButton.setGeometry(QtCore.QRect(1030, 40, 112, 32))
		self.projectInfoButton.setObjectName("projectInfoButton")
		self.label_3 = QtWidgets.QLabel(Form)
		self.label_3.setGeometry(QtCore.QRect(780, 20, 281, 20))
		self.label_3.setObjectName("label_3")
		self.projectInfoEntry = QtWidgets.QLineEdit(Form)
		self.projectInfoEntry.setGeometry(QtCore.QRect(780, 40, 241, 32))
		self.projectInfoEntry.setObjectName("projectInfoEntry")
		self.sampleInfoButton = QtWidgets.QPushButton(Form)
		self.sampleInfoButton.setGeometry(QtCore.QRect(1030, 110, 112, 32))
		self.sampleInfoButton.setObjectName("sampleInfoButton")
		self.label_4 = QtWidgets.QLabel(Form)
		self.label_4.setGeometry(QtCore.QRect(780, 90, 281, 20))
		self.label_4.setObjectName("label_4")
		self.sampleInfoEntry = QtWidgets.QLineEdit(Form)
		self.sampleInfoEntry.setGeometry(QtCore.QRect(780, 110, 241, 32))
		self.sampleInfoEntry.setObjectName("sampleInfoEntry")
		self.readsInfoButton = QtWidgets.QPushButton(Form)
		self.readsInfoButton.setGeometry(QtCore.QRect(1030, 190, 112, 32))
		self.readsInfoButton.setObjectName("readsInfoButton")
		self.label_5 = QtWidgets.QLabel(Form)
		self.label_5.setGeometry(QtCore.QRect(780, 170, 281, 20))
		self.label_5.setObjectName("label_5")
		self.readsInfoEntry = QtWidgets.QLineEdit(Form)
		self.readsInfoEntry.setGeometry(QtCore.QRect(780, 190, 241, 32))
		self.readsInfoEntry.setObjectName("readsInfoEntry")
		self.label_6 = QtWidgets.QLabel(Form)
		self.label_6.setGeometry(QtCore.QRect(470, 280, 141, 20))
		self.label_6.setObjectName("label_6")
		self.usernameEntry = QtWidgets.QLineEdit(Form)
		self.usernameEntry.setGeometry(QtCore.QRect(470, 300, 171, 32))
		self.usernameEntry.setObjectName("usernameEntry")
		self.passwordEntry = QtWidgets.QLineEdit(Form)
		self.passwordEntry.setGeometry(QtCore.QRect(660, 300, 171, 32))
		self.passwordEntry.setEchoMode(QtWidgets.QLineEdit.Password)
		self.passwordEntry.setObjectName("passwordEntry")
		self.label_7 = QtWidgets.QLabel(Form)
		self.label_7.setGeometry(QtCore.QRect(660, 280, 141, 20))
		self.label_7.setObjectName("label_7")
		self.newProjectCombo = QtWidgets.QComboBox(Form)
		self.newProjectCombo.setGeometry(QtCore.QRect(530, 40, 90, 28))
		self.newProjectCombo.setObjectName("newProjectCombo")
		self.newProjectCombo.addItem("")
		self.newProjectCombo.addItem("")
		self.label_8 = QtWidgets.QLabel(Form)
		self.label_8.setGeometry(QtCore.QRect(530, 20, 161, 20))
		self.label_8.setObjectName("label_8")
		self.newSampleCombo = QtWidgets.QComboBox(Form)
		self.newSampleCombo.setGeometry(QtCore.QRect(530, 110, 90, 28))
		self.newSampleCombo.setObjectName("newSampleCombo")
		self.newSampleCombo.addItem("")
		self.newSampleCombo.addItem("")
		self.label_9 = QtWidgets.QLabel(Form)
		self.label_9.setGeometry(QtCore.QRect(530, 90, 181, 20))
		self.label_9.setObjectName("label_9")
		self.testSubmissionCombo = QtWidgets.QComboBox(Form)
		self.testSubmissionCombo.setGeometry(QtCore.QRect(530, 190, 90, 28))
		self.testSubmissionCombo.setObjectName("testSubmissionCombo")
		self.testSubmissionCombo.addItem("")
		self.testSubmissionCombo.addItem("")
		self.label_10 = QtWidgets.QLabel(Form)
		self.label_10.setGeometry(QtCore.QRect(530, 170, 161, 20))
		self.label_10.setObjectName("label_10")
		self.label_11 = QtWidgets.QLabel(Form)
		self.label_11.setGeometry(QtCore.QRect(20, 340, 76, 20))
		self.label_11.setObjectName("label_11")
		self.logArea = QtWidgets.QTextEdit(Form)
		self.logArea.setGeometry(QtCore.QRect(20, 360, 911, 241))
		self.logArea.setObjectName("logArea")
		self.label_12 = QtWidgets.QLabel(Form)
		self.label_12.setGeometry(QtCore.QRect(950, 350, 211, 221))
		self.label_12.setText("")
		self.label_12.setPixmap(QtGui.QPixmap(installationDirectory+"src/GUI/IconsFinal/DB.jpg"))
		self.label_12.setObjectName("label_12")
		self.pushButton_6 = QtWidgets.QPushButton(Form)
		self.pushButton_6.setGeometry(QtCore.QRect(950, 570, 201, 32))
		self.pushButton_6.setObjectName("pushButton_6")
		self.selectedFilesArea = QtWidgets.QTextEdit(Form)
		self.selectedFilesArea.setGeometry(QtCore.QRect(20, 40, 301, 291))
		self.selectedFilesArea.setObjectName("selectedFilesArea")

		self.retranslateUi(Form)
		QtCore.QMetaObject.connectSlotsByName(Form)

		#Addition to GUI generated by Qt designer
		self.projectInfoButton.clicked.connect(self.selectProjectTable)
		self.sampleInfoButton.clicked.connect(self.selectSampleTable)
		self.readsInfoButton.clicked.connect(self.selectReadsInfo)
		self.pushButton_6.clicked.connect(self.selectFiles)


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
		if str(self.selectedFilesArea.toPlainText) == "":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("No fastq files has been selected.")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should select at least one paired end reads dataset to submit to ENA\n ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return
		
		projInfo = self.projectInfoEntry.text()
		if projInfo == "" and self.newProjectCombo.currentText()=="Yes":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid project info table should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("It seems like you chose to create a new project (aka study) on the ENA database. You should provide a table reporting the project information to do so (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		samInfo = self.sampleInfoEntry.text()
		if samInfo == "" and self.newSampleCombo.currentText()=="Yes":
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid sample info table should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("It seems like you chose to create a new sample for each dataset on the ENA database. You should provide a table reporting the samples information to do so (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		reInfo = self.readsInfoEntry.text()
		if reInfo == "" :
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid reads info table should be provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should provide a table with all the required read datasets information in order to submit them to ENA (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		usN = self.usernameEntry.text()
		if usN == "" :
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid ENA username was not provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should insert your ENA username and password (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		usN = self.passwordEntry.text()
		if usN == "" :
			msg = QMessageBox()
			msg.setIcon(QMessageBox.Warning)
			msg.setText("A valid ENA password was not provided")
			msg.setWindowTitle("Warning")
			msg.setDetailedText("You should insert your ENA username and password (see manual) ")
			msg.setStandardButtons(QMessageBox.Ok)
			msg.exec_()
			return

		#*******************************************************************
		#********************* Main algorithm reads submission  ************
		#*******************************************************************
		isTest = self.testSubmissionCombo.currentText()

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

		for dataset in self.onlyfiles:

			sampleName = self.getPrefix((dataset[0].split("/"))[-1])
			self.refreshTextArea(sampleName)
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

		for item in self.onlyfiles:
			item2 = self.getPrefix((item[0].split("/"))[-1])
			#print(item,sampleAccession[item],runAccession[item],experimentAccession[item])
			outfile.write(item2+"\t"+sampleAccession[item2]+"\t"+runAccession[item2]+"\t"+experimentAccession[item2]+"\n")

		outfile.close()


	def retranslateUi(self, Form):
		_translate = QtCore.QCoreApplication.translate
		Form.setWindowTitle(_translate("Form", "Sequencing data submission tool"))
		self.label.setText(_translate("Form", "Samples to submit"))
		self.sampleToSubmitButton.setText(_translate("Form", "Open file"))
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

