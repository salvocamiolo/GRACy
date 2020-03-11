#!/usr/bin/python
installationDirectory = "/home/gracy/Desktop/GRACy/"


try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

from PIL import ImageTk, Image
import sys
import tkFont
import os



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

        def launchQC():
            os.system("python "+installationDirectory+"qualityCheck/readsQualityCheck.py "+installationDirectory+" &")

        def denovoAssembly():
            os.system("python "+installationDirectory+"assembly/hcmvAssembly.py "+installationDirectory+" &")

        def genotyping():
            os.system("python "+installationDirectory+"genotyping/hcmvGenotypingGUI.py "+installationDirectory+" &")

        def annotation():
            os.system("python "+installationDirectory+"annotation/hcmvAnnonation.py "+installationDirectory+" &")

        def snpAnalysis():
            os.system("python "+installationDirectory+"snpAnalysis/hcmvSNPAnalysisGUI.py "+installationDirectory+" &")

        def dbSubmission():
            os.system("python "+installationDirectory+"databaseSubmission/hcmvDataSubmission.py "+installationDirectory+" &")




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

        w=1000
        h=700
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y))
        top.title("GRACy")
        top.configure(highlightcolor="black",background="white")

        #image = ImageTk.PhotoImage(file=installationDirectory+"resources/Medicon_virus.jpg")

        #gmail=ImageTk.PhotoImage(file=installationDirectory+'resources/Medicon_virus_times.jpg')
        #self.lab=tk.Label(image=gmail)
        #self.lab.photo=gmail
        #self.lab.pack()

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/filtering.jpg")
        photo1 =ImageTk.PhotoImage(image)
        readsQCButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=launchQC,image=photo1)
        readsQCButton.place(x=30,y=30)
        readsQCButton.image = photo1

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/Classification.jpg")
        photo1 =ImageTk.PhotoImage(image)
        gentypingButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=genotyping,image=photo1)
        gentypingButton.place(x=280,y=30)
        gentypingButton.image = photo1
	
	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/Assembly.jpg")
        photo1 =ImageTk.PhotoImage(image)
        assemblyButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=denovoAssembly,image=photo1)
        assemblyButton.place(x=530,y=30)
        assemblyButton.image = photo1

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/Annotation.jpg")
        photo1 =ImageTk.PhotoImage(image)
        annotationButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=annotation,image=photo1)
        annotationButton.place(x=780,y=30)
        annotationButton.image = photo1

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/SNPcalling.jpg")
        photo1 =ImageTk.PhotoImage(image)
        snpAnalysisButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=snpAnalysis,image=photo1)
        snpAnalysisButton.place(x=30,y=280)
        snpAnalysisButton.image = photo1

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/DB.jpg")
        photo1 =ImageTk.PhotoImage(image)
        dbSubmissionButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=dbSubmission,image=photo1)
        dbSubmissionButton.place(x=280,y=280)
        dbSubmissionButton.image = photo1

	image = Image.open(installationDirectory+"resources/imagesInterface/IconsFinal/IconGRACy.jpg")
        photo1 =ImageTk.PhotoImage(image)
        logoLabel = tk.Label(top, compound=tk.TOP,height=150,width=150,image=photo1,borderwidth=0,highlightthickness=0)
        logoLabel.place(x=820,y=520)
        logoLabel.image = photo1
	
        #readsQCButton = tk.Button(top, command=launchQC)
        #readsQCButton.place(x=535,y=30,height=30,width=180)
        #readsQCButton.configure(text= "Reads QC",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')

        #assemblyButton = tk.Button(top,command=denovoAssembly)
        #assemblyButton.place(x=640,y=265,height=30,width=180)
        #assemblyButton.configure(text= "De novo assembly",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')

        #gentypingButton = tk.Button(top,command=genotyping)
        #gentypingButton.place(x=600,y=412,height=30,width=150)
        #gentypingButton.configure(text= "Genotyping",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')

        #annotationButton = tk.Button(top,command = annotation)
        #annotationButton.place(x=200,y=495,height=30,width=150)
        #annotationButton.configure(text= "Annotation",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')

        #dbSubmissionButton = tk.Button(top,command = dbSubmission)
        #dbSubmissionButton.place(x=65,y=167,height=30,width=150)
        #dbSubmissionButton.configure(text= "DB submission",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')

        #snpAnalysisButton = tk.Button(top,command=snpAnalysis)
        #snpAnalysisButton.place(x=65,y=358,height=30,width=150)
        #snpAnalysisButton.configure(text= "SNP analysis",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98",activebackground='#568F98')




if __name__ == '__main__':
    vp_start_gui()
