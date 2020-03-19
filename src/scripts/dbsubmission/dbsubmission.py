import tkinter as tk
from  tkinter import  ttk
import tkinter
from tkinter import filedialog as tkFileDialog
from PIL import ImageTk, Image
import sys
from tkinter import messagebox

installationDirectory = sys.argv[1]
import os
import time

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
    
    top = Toplevel1 (w)
    
    return (w, top)

def destroy_Toplevel1():
    global w
    w.destroy()
    w = None

class Toplevel1:
    def __init__(self, top=None):


        sampleAccession = {}
        runAccession = {}
        projectAccessionNumber = ""
        experimentAccession = {}
        sampleList = []
        
        def exitProgram():
            exit()

        def openRunToSubmit():
            toSubmit = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
            self.runToSubmitFileEntry.delete(0,tk.END)
            self.runToSubmitFileEntry.insert(0,toSubmit)
        
        def openFastqFolder():
            fastqFodler = tkFileDialog.askdirectory(initialdir = "./",title = "Select folder")
            self.fastqFolderEntry.delete(0,tk.END)
            self.fastqFolderEntry.insert(0,fastqFodler)

        def openFastqInfoFile():
            fastqInfo = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
            self.fastqInfoEntry.delete(0,tk.END)
            self.fastqInfoEntry.insert(0,fastqInfo)

        def openProjectInfoFile():
            projectInfo = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
            self.projectInfoEntry.delete(0,tk.END)
            self.projectInfoEntry.insert(0,projectInfo)

        def openSampleInfoFile():
            sampleInfo = tkFileDialog.askopenfilename(initialdir = "./",title = "Select file")
            self.sampleInfoEntry.delete(0,tk.END)
            self.sampleInfoEntry.insert(0,sampleInfo)

        def createNewProject(projectFile):
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
                if self.testSubmissionchkValue.get()==True:
                    os.system("curl -u "+self.usernameEntry.get()+":"+self.passwordEntry.get()+" -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt")
                else:
                    MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
                    if MsgBox == 'yes':
                        os.system("curl -u "+self.usernameEntry.get()+":"+self.passwordEntry.get()+" -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt")
                
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

        def createNewSample(sampleName):
            if not sampleName in sampleAccession:
                sampleAccession[sampleName] = ""
            infile = open(self.sampleInfoEntry.get())

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
                    if self.testSubmissionchkValue.get()==True:
                        os.system("curl -u "+self.usernameEntry.get()+":"+self.passwordEntry.get()+" -F \"SUBMISSION=@submission.xml\" -F \"SAMPLE=@"+alias+"_sample.xml\" \"https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/\" >sampleReceipt.txt")
                    else:
                        MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
                        if MsgBox == 'yes':
                            os.system("curl -u "+self.usernameEntry.get()+":"+self.passwordEntry.get()+" -F \"SUBMISSION=@submission.xml\" -F \"SAMPLE=@"+alias+"_sample.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >sampleReceipt.txt")

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

        


        
        
        
        def readsSubmission():
            if self.runToSubmitFileEntry.get()=="Please select a file...." or self.fastqFolderEntry.get()=="Please select a folder...." or self.fastqInfoEntry.get()=="Please select a file...." or self.usernameEntry.get()=="myENA_Username" or self.passwordEntry.get()=="":
                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "The following fields must be filled:\n")
                self.logArea.insert(tk.END, "Samples to submit\nFastq folder\nFastq info folder\nENA username\nPassword\n\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()
            else:
                #*******************************************************************
                #********************* Main algorithm reads submission  ************
                #*******************************************************************
                isTest = self.testSubmissionchkValue.get() 
                inputFolder = self.fastqFolderEntry.get()
                filesToSubmit = self.runToSubmitFileEntry.get()
                createProject = self.createProjectchkValue.get()
                createSample = self.createSampleschkValue.get()
                projectInfo = self.projectInfoEntry.get()
                sampleInfo = self.sampleInfoEntry.get()
                readsInfo = self.fastqInfoEntry.get()
                os.system("cp "+installationDirectory+"src/scripts/dbsubmission/utils/submission.xml .")

                

                self.logArea.configure(state='normal')
                self.logArea.insert(tk.END, "isTest: "+str(isTest)+"\n"+"inputFolder: "+str(inputFolder)+"\nfilesToSubmit: "+str(filesToSubmit)+"\nCreatePoject: "+str(createProject)+"\nCreateSample: "+str(createSample)+"\nProjectInfo: "+str(projectInfo)+"\nSampleInfo: "+str(sampleInfo)+"\nReadsInfo: "+str(readsInfo)+"\n\n")
                self.logArea.see(tk.END)
                self.logArea.configure(state='disabled')
                self.logArea.update()

                if createProject==True:
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "Creating a new Project....\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()
                    projectAccessionNumber = createNewProject(projectInfo)
                    self.logArea.configure(state='normal')
                    self.logArea.insert(tk.END, "A project was created on ENA with the accession number: "+projectAccessionNumber+"\n\n")
                    self.logArea.see(tk.END)
                    self.logArea.configure(state='disabled')
                    self.logArea.update()

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
                    if createSample == True:
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Creating a new sample for dataset "+sampleName+"....\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()
                        sampleAccessionNumber = createNewSample(sampleName)
                        self.logArea.configure(state='normal')
                        self.logArea.insert(tk.END, "Sample created for dataset "+sampleName+" with accession number: "+sampleAccessionNumber+"\n")
                        self.logArea.see(tk.END)
                        self.logArea.configure(state='disabled')
                        self.logArea.update()

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

                        if createProject == True:
                            projectAccessionField = projectAccessionNumber
                        if createSample ==True:
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
                            
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Submitting reads for sample "+sampleName+"....\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()

                            if self.testSubmissionchkValue.get()==True:
                                os.system(installationDirectory+"src/conda/bin/java -Xmx2048m -jar "+installationDirectory+"src/scripts/dbsubmission/utils/webin-cli-2.2.0.jar  -context reads -userName "+self.usernameEntry.get()+" -password "+self.passwordEntry.get()+"  -manifest "+sampleName+"_manifestFile.txt -test -submit >fastqReceipt")
                            else:
                                MsgBox = tk.messagebox.askquestion("You are going to submit your data to ENA. This is not a test submission. Do you with to continue?")
                                if MsgBox == 'yes':
                                    os.system(installationDirectory+"src/conda/bin/java -Xmx2048m -jar "+installationDirectory+"src/scripts/dbsubmission/utils/webin-cli-2.2.0.jar -context reads -userName "+self.usernameEntry.get()+" -password "+self.passwordEntry.get()+"  -manifest "+sampleName+"_manifestFile.txt -submit >fastqReceipt")
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
                            self.logArea.configure(state='normal')
                            self.logArea.insert(tk.END, "Fastq files submitted!\n")
                            self.logArea.insert(tk.END, "Run accession number: "+runAccession[sampleName]+"\n")
                            self.logArea.insert(tk.END, "Experiment accession number: "+experimentAccession[sampleName]+"\n")
                            self.logArea.see(tk.END)
                            self.logArea.configure(state='disabled')
                            self.logArea.update()
                            
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
        h=550
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("Data submission tool")
        top.configure(highlightcolor="black",background="white")

        self.runToSubmitFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.runToSubmitFileLabel.configure(text="Samples to submit")
        self.runToSubmitFileLabel.place(x=20,y=20,width=130,height=20)

        self.runToSubmitFileEntry = tk.Entry(top)
        self.runToSubmitFileEntry.place(x=20,y=40,width=300,height=30)
        self.runToSubmitFileEntry.insert(0,"Please select a file....")

        self.runToSubmitButton = tk.Button(command=openRunToSubmit,background="#204949",foreground="white")
        self.runToSubmitButton.place(x=330,y=40,width=100,height=30)
        self.runToSubmitButton.configure(text="Open file")

        self.fastqFolderLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.fastqFolderLabel.configure(text="Fastq folder")
        self.fastqFolderLabel.place(x=20,y=80,width=80,height=20)

        self.fastqFolderEntry = tk.Entry(top)
        self.fastqFolderEntry.place(x=20,y=100,width=300,height=30)
        self.fastqFolderEntry.insert(0,"Please select a folder....")

        self.fastqFolderButton = tk.Button(top,command=openFastqFolder,background="#204949",foreground="white")
        self.fastqFolderButton.place(x=330,y=100,width=100,height=30)
        self.fastqFolderButton.configure(text="Open folder")

        self.fastqInfoLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.fastqInfoLabel.configure(text="Fastq info table")
        self.fastqInfoLabel.place(x=20,y=140,width=110,height=20)

        self.fastqInfoEntry = tk.Entry(top)
        self.fastqInfoEntry.place(x=20,y=160,width=300,height=30)
        self.fastqInfoEntry.insert(0,"Please select a file....")

        self.fastqInfoFileButton = tk.Button(top,command=openFastqInfoFile,background="#204949",foreground="white")
        self.fastqInfoFileButton.place(x=330,y=160,width=100,height=30)
        self.fastqInfoFileButton.configure(text="Open file")


        self.projectInfoLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.projectInfoLabel.configure(text="Project info table")
        self.projectInfoLabel.place(x=470,y=20,width=110,height=20)

        self.projectInfoEntry = tk.Entry(top)
        self.projectInfoEntry.place(x=470,y=40,width=300,height=30)
        self.projectInfoEntry.insert(0,"Please select a file....")

        self.projectInfoButton = tk.Button(top,command=openProjectInfoFile,background="#204949",foreground="white")
        self.projectInfoButton.place(x=780,y=40,width=100,height=30)
        self.projectInfoButton.configure(text="Open file")


        self.sampleInfoLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.sampleInfoLabel.configure(text="Samples info Table")
        self.sampleInfoLabel.place(x=470,y=80,width=130,height=20)

        self.sampleInfoEntry = tk.Entry(top)
        self.sampleInfoEntry.place(x=470,y=100,width=300,height=30)
        self.sampleInfoEntry.insert(0,"Please select a file....")

        self.sampleInfoButton = tk.Button(top,command=openSampleInfoFile,background="#204949",foreground="white")
        self.sampleInfoButton.place(x=780,y=100,width=100,height=30)
        self.sampleInfoButton.configure(text="Open file")


        self.usernameLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.usernameLabel.configure(text="Username")
        self.usernameLabel.place(x=470,y=140,width=75,height=20)

        self.usernameEntry = tk.Entry(top,justify='right')
        self.usernameEntry.place(x=470,y=160,width=200,height=30)
        self.usernameEntry.insert(0,"myENA_Username")

        self.passwordLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.passwordLabel.configure(text="Password")
        self.passwordLabel.place(x=680,y=140,width=75,height=20)

        self.passwordEntry = tk.Entry(top,justify='right', show="*")
        self.passwordEntry.place(x=680,y=160,width=200,height=30)
        self.passwordEntry.insert(0,"")

        self.createProjectchkValue = tk.BooleanVar() 
        self.createProjectchkValue.set(True)
        self.createProjectCheckButton = tk.Checkbutton(top,variable=self.createProjectchkValue,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.createProjectCheckButton.place(x=20,y=230,height=20,width=150)
        self.createProjectCheckButton.configure(text="Create new Project")
        #createProjectCheckButton.bind('<Button-1>',selectKraken)

        self.createSampleschkValue = tk.BooleanVar() 
        self.createSampleschkValue.set(True)
        self.createSamplesCheckButton = tk.Checkbutton(top,variable=self.createProjectchkValue,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.createSamplesCheckButton.place(x=200,y=230,height=20,width=170)
        self.createSamplesCheckButton.configure(text="Create new Sample(s)")

        self.testSubmissionchkValue = tk.BooleanVar() 
        self.testSubmissionchkValue.set(True)
        self.testSubmissionCheckButton = tk.Checkbutton(top,variable=self.testSubmissionchkValue,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.testSubmissionCheckButton.place(x=380,y=230,height=20,width=150)
        self.testSubmissionCheckButton.configure(text="Test submission")


        self.runButton = tk.Button(top,command=readsSubmission,background="#204949",foreground="white")
        self.runButton.place(x=680,y=500,width=200,height=30)
        self.runButton.configure(text="Submit")

        #self.exitButton = tk.Button(top,command=exitProgram)
        #self.exitButton.place(x=660,y=220,width=100,height=30)
        #self.exitButton.configure(text="Exit")
       

        self.logFileLabel = tk.Label(top,highlightthickness=0,background="white",borderwidth=0,foreground="#204949")
        self.logFileLabel.place(x=20,y=270,height=20,width=100)
        self.logFileLabel.configure(text="Log window")

        self.logFrame = tk.Frame(top)
        self.logFrame.place(x=20, y=290, height=240, width=640)
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(borderwidth="2")
        self.logFrame.configure(relief='groove')
        self.logFrame.configure(width=125)
        self.logArea = tk.Text(top,state='disabled')
        self.logArea.place(x=25,y=295,height=230, width=630)
        self.logArea.configure(background="white",borderwidth=5)
        self.logArea.configure(selectbackground="#c4c4c4")

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/DB.jpg")
        photo1 =ImageTk.PhotoImage(image)
        logoLabel = tk.Label(top, compound=tk.TOP,height=200,width=200,image=photo1,borderwidth=0,highlightthickness=0)
        logoLabel.place(x=680,y=290)
        logoLabel.image = photo1


        

if __name__ == '__main__':
    vp_start_gui()
