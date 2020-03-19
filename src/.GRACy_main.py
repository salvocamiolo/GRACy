

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
from tkinter import font as tkFont
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
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/readsFiltering/readsFilteringQt.py "+installationDirectory+" &")

		def denovoAssembly():
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/assembly/assemblyQt.py "+installationDirectory+" &")

		def genotyping():
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/genotyping/genotypingQt.py "+installationDirectory+" &")

		def annotation():
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/annotation/annotationQt.py "+installationDirectory+" &")

		def snpAnalysis():
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/snpCalling/snpCallingQt.py "+installationDirectory+" &")

		def dbSubmission():
			os.system(installationDirectory+"src/conda/bin/python "+installationDirectory+"src/scripts/dbsubmission/dbsubmissionQt.py "+installationDirectory+" &")




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

		x = (ws/2) - (w/2)
		y = (hs/2) - (h/2)

		top.geometry('%dx%d+%d+%d' % (w, h, x, y))
		top.title("GRACy")
		top.configure(highlightcolor="black",background="white")

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/filtering.jpg")
		photo1 =ImageTk.PhotoImage(image)
		readsQCButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=launchQC,image=photo1)
		readsQCButton.place(x=30,y=30)
		readsQCButton.image = photo1

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Classification.jpg")
		photo1 =ImageTk.PhotoImage(image)
		gentypingButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=genotyping,image=photo1)
		gentypingButton.place(x=280,y=30)
		gentypingButton.image = photo1
	
		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Assembly.jpg")
		photo1 =ImageTk.PhotoImage(image)
		assemblyButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=denovoAssembly,image=photo1)
		assemblyButton.place(x=530,y=30)
		assemblyButton.image = photo1

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/Annotation.jpg")
		photo1 =ImageTk.PhotoImage(image)
		annotationButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=annotation,image=photo1)
		annotationButton.place(x=780,y=30)
		annotationButton.image = photo1

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/SNPcalling.jpg")
		photo1 =ImageTk.PhotoImage(image)
		snpAnalysisButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=snpAnalysis,image=photo1)
		snpAnalysisButton.place(x=30,y=280)
		snpAnalysisButton.image = photo1

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/DB.jpg")
		photo1 =ImageTk.PhotoImage(image)
		dbSubmissionButton = tk.Button(top, compound=tk.TOP,height=200,width=200,command=dbSubmission,image=photo1)
		dbSubmissionButton.place(x=280,y=280)
		dbSubmissionButton.image = photo1

		image = Image.open(installationDirectory+"src/GUI/IconsFinal/IconGRACy.jpg")
		photo1 =ImageTk.PhotoImage(image)
		logoLabel = tk.Label(top, compound=tk.TOP,height=150,width=150,image=photo1,borderwidth=0,highlightthickness=0)
		logoLabel.place(x=820,y=520)
		logoLabel.image = photo1
	



if __name__ == '__main__':
	vp_start_gui()
