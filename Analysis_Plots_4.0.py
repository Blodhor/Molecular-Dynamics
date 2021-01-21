import numpy as np
import matplotlib.pyplot as plt

class Analysis_plot:
	def __init__(self, type = 'one', names=[('file.dat','analysis_title')], analysisType = 'rmsd', largerYaxis = False, 
	Yenlarger = 2, frameToTime=False, frameStep = 5*10**4, timeStep = 0.004, nanosec = False, suptitle='Titulo geral', 
	labelpx = 35.0, labelpy = 0.50, dpi = 100):
		'''Parameters:
		
		type: Plot 'one' dataset, 'two' datasets beside each other or 'four' datasets in a matrix fashion.

		names: Data Files, the format expected is .dat, if any other is used do not expect a pretty result.
		
		analysisType: Y-axis label between: rmsd; rmsf; radgyr.
		
		largerYaxis: Boolean argument. Sets if the Y-axis  will enlarge or not.
		
		Yenlarger: Sets the Y-axis max value, mutiplying itself to max(Y-dataset). 
		
		frameToTime: Boolean argument. For rmsd and radgyr analysis. If True X-axis will change from frame to time (pico seconds).
			
		frameStep: Integer value. Amount of 'timeStep's used in the production. Default 50 000 (gpu-high mode from ASM.py).
		
		timeStep: 0.002 ps or 0.004 ps for HMR. Default 0.004 ps.

		nanosec: Boolean argument. If True "time" marks on X-axis will be on nanoseconds.
		
		suptitle: Title of the whole plot. Giving titles to every subplot when type is 'four' decreases the pict. quality.
		
		labelpx: X position of internal labelling.
		
		labelpy: Y position of internal labelling.
		
		dpi: Sets the picture quality.'''

		self.X = []
		self.Y = []
		self.frameToTime = frameToTime
		self.frameStep = frameStep
		self.timeStep = timeStep
		self.nanosec = nanosec
		self.dpi = dpi
		self.suptitle = suptitle
		self.labelpx = labelpx
		self.labelpy = labelpy
		 
		if 'radgyr' in analysisType:
			self.ana_type = 'Raio de giro (Angstrom)'
		else:
			self.ana_type = analysisType.upper()

		XTlabels = self.XY(files=names)

		if XTlabels != -1 and type == 'one':

			self.plot_one(X=self.X[0], Y=self.Y[0], Xaxis=XTlabels[0][0], 
			name= XTlabels[0][1], EnlargeYaxis=largerYaxis, Ymax=Yenlarger)

		elif XTlabels != -1 and type == 'two':

			self.plot_two(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[0][1], XTlabels[1][1]])

		elif XTlabels != -1 and type == 'four':

			self.plot_four(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[0][1], XTlabels[1][1], XTlabels[2][1], XTlabels[3][1]])

		else:
			print("Type must be one, two or four!\n")

	def XY(self, files = [('file1.dat','title1'), ('file2.dat','title2')]):
		''' Gets X and Y datasets and returns its respectives x-label and title as a list of tuples.'''

		if len(files) == 1 or len(files) == 2 or len(files) == 4:

			xname = 'T'
			XTlabel = []

			for fs in files:
				f = open(fs[0])
				t = f.readlines()
				f.close()

				plot_title = fs[1]

				x = []
				y = []

				for i in t:
					data = i.split()
					if i != '' and i[0] == '#':
						if 'RMSF' in self.ana_type:
							xname = 'Resíduo'
						elif not self.frameToTime:
							xname = data[0][1:]
							if 'frame' in xname.lower():
								xname = 'Instância'
						else:
							if self.nanosec:
								xname = 'Tempo (ns)'
							else:
								xname = 'Tempo (ps)'
					elif len(data) == 2:
						if 'RMSF' in self.ana_type:
							x.append( int(float(data[0])) )
						elif self.frameToTime:
							if self.nanosec:
								x.append( float(data[0])* self.frameStep * self.timeStep/1000)
							else:
								x.append( float(data[0])* self.frameStep * self.timeStep)
						else:
							x.append( float(data[0]) )

						y.append( float(data[1]) )

				self.X.append(x)
				self.Y.append(y)
				XTlabel.append( (xname,plot_title) )

			return XTlabel
		else:
			return -1

	def plot_one(self, X = [], Y = [], Xaxis = "Frames", name = "Título", EnlargeYaxis = False, Ymax = 2):
		'''Plots one X-Y dataset.'''
		plt.figure(dpi=self.dpi)

		if EnlargeYaxis:
			#lines below set x and y axis limits
			axes = plt.axes()
			axes.set_xlim([0,X[len(X)-1]])
			axes.set_ylim([0,Ymax*Y[len(Y)-1]])
			#if for some reason the plot differente than expected you can try changing the x,y spacements.
			#axes.set_xticks(np.arange(0,X[len(X)-1],X[len(X)-1]/8.0)) #numpy bugged hard with the last windows patch
			#axes.set_yticks(np.arange(0,Ymax*Y[len(Y)-1],Ymax*Y[len(Y)-1]/5.0))

		plt.plot(X,Y)
		#It's possible to make a log, dilog,... 
		#plt.xscale('linear') #'log'
		#plt.yscale('linear')

		plt.title(name)
		plt.ylabel(self.ana_type)
		plt.xlabel(Xaxis)
		plt.show()

	def plot_two(self, X = [[], []], Y = [[], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''
		
		fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, dpi=self.dpi)

		ax1.plot(X[0],Y[0])
		ax1.set_title(name[0])

		ax2.plot(X[1],Y[1])
		ax2.set_title(name[1])

		for ts in [ax1,ax2]:
			ts.set(xlabel=Xaxis, ylabel=self.ana_type)
		
		#pyplot.subplots can hide redundant axes
		for ts in [ax1,ax2]:
			ts.label_outer()
		
		plt.show()

	def plot_four(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''
		
		plt.figure(dpi=self.dpi)

		ax1 = plt.subplot(221)
		plt.plot(X[0],Y[0])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[0], color='darkblue')
		plt.ylabel(self.ana_type)
		#plt.xlabel(Xaxis) #Showing this will interfere with the quality of the pict.
		plt.grid(True)
		
		ax2 = plt.subplot(222, sharey=ax1)		
		plt.plot(X[1],Y[1])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[1], color='darkblue')
		#plt.ylabel(self.ana_type) #Showing this will interfere with the quality of the pict.
		#plt.xlabel(Xaxis) #Showing this will interfere with the quality of the pict.
		plt.grid(True)

		ax3 = plt.subplot(223, sharex=ax1)
		plt.plot(X[2],Y[2])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[2], color='darkblue')
		plt.ylabel(self.ana_type)
		plt.xlabel(Xaxis)
		plt.grid(True)
		
		ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
		plt.plot(X[3],Y[3])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[3], color='darkblue')
		plt.xlabel(Xaxis)
		plt.grid(True)
		
		plt.suptitle(self.suptitle)
		plt.show()

if __name__ == "__main__":
	import sys
	arg = sys.argv
	version_n = 4.0
	version = ''' Analysis plot - Pyplot extension for RMSD, RMSF and Radgyr data analysis.

Python3 modules needed: 
	numpy; 
	matplotlib.

Data file structure expected:

----------------------------------------------------------------------------------
| ph10-rmsd.dat (200 frames) | ph10-radgyr.dat       | ph10-rmsf.dat (res 6-260) |
-----------------------------|-----------------------|----------------------------
| #Frame          reference  | #Frame    ph10-radgyr | #Res        AtomicFlx     |
|       1       0.0000       |        1      16.4009 |    6.000       1.0428     |
|       2       0.6861       |        2      16.4335 |    7.000       0.8080     |
|       3       0.6666       |        3      16.4809 | .......                   |
|                                   ......                                       |
|       199       0.8485     |      199      16.5353 | 259.000       0.4021      |
|       200       0.9146     |      200      16.5919 | 260.000       0.4458      |
----------------------------------------------------------------------------------



Autor: Braga, B. C.
e-mail: bruno.braga@ufms.br '''

	#default keys
	version_only = False
	inst_only    = False
	tpe          = 'one'
	anatp        = 'rmsd'
	File         = [('Gpu-Ultra_cphmd-petase-water/ph7.00_rmsd.dat','pH=7.00')]
	supertitle   = 'Petase nativa - Produção'
	labx, laby   = (30.5, 0.51) # rmsf: (150, 2.25) # radgyr: (160, 16.61)
	fram2time    = False
	framstp      = 25*10**4 # Ultra
	nano         = False
	dpi          = 100

	flags = ["&","-v","--version","-h","--help",'-type', '-i', '-anatp','-stitle','-fram2time', '-framstp','-nanosec','-dpi','-lblcrd']

	cut =0 # counter for input flags

	for i in range(len(arg)):
		if cut == i:
			if arg[i].lower() == "&":
				break
			elif arg[i].lower() == "-v" or arg[i].lower() == "--version":
				version_only = True
				print("Current version: %s\n"%version_n)
				print(version)
				break
			elif arg[i].lower() == "-h" or arg[i].lower() == "--help":
				inst_only = True
				print("Welcome to Analysis plot %s:\n"%version_n)
				print("\nUsage:\n\t-v or --version\t\tprints current version, the python libraries needed and the data format expected.\n")
				print("\t-h or --help\t\tprints this message.\n")
				print("\t-type\tQuantity of files used, for data comparison.\n\t\tone: Normal plot with one data file.\n\t\ttwo: Plots two data files with the same axis information (ex:RMSD for pH7 and pH8).\n\t\tfour: Plots four data files with the same axis information (ex:RMSD for pH5, pH6, pH7, pH8).\n")
				print("\t-anatp\t\tAnalysis type:\n\t\trmsd;\n\t\trmsf;\n\t\tradgyr.\n")
				print("\t-i\t\tinput data file(s) with a name for the plot (separated by space).\n\t\t\tEx: -i ph7.00_rmsd.dat pH=7.00\n")
				print("\t-stitle\t\t(Valid only for type four) Title for comparison plot.\n")
				print("\t-lblcrd\t(Valid only for type four) Label coords.\n")
				print("\t-fram2time\t(For anatp = rmsd or radgyr) Sets X-axis will change from frame to time (pico seconds). Must inform frame-step conversion.\n\t\t\tEx:-fram2time 10000\n")
				print("\t-nanosec\t\t(For anatp = rmsd or radgyr) Sets time intervals on X-axis to nanoseconds.\n")
				print("\t-dpi\t\tSets the plot quality (recommended to use only with type=one). Default: 100\n")
				print("\nExamples:\n\t$ python3 Analysis_Plots_4.0.py -type one -anatp rmsd -i ph7.00_rmsd.dat Production pH=7.00 -fram2time 10000 -dpi 120\n")
				print("\n\t$ python3 Analysis_Plots_4.0.py -type two -anatp rmsf -i ph7.00_rmsf.dat Production pH=7.00 ph8.00_rmsf.dat Production pH=8.00\n")
				print("\n\t$ python3 Analysis_Plots_4.0.py -type four -stitle Production -anatp radgyr -lblcrd (16.61,200) -i ph7.00_radgyr.dat pH=7.00 ph8.00_radgyr.dat pH=8.00 ph9.00_radgyr.dat pH=9.00 ph10.00_radgyr.dat pH=10.00 -fram2time 10000 -nanosec\n")
				break

			elif arg[i].lower() == "-type":
				tpe = arg[i+1]
				continue
			elif arg[i].lower() == "-anatp":
				anatp = arg[i+1]
				continue
			elif arg[i].lower() == "-i":
				File = []
				#File.append( (arg[i+1], arg[i+2]) ) # arg[i+2] != 'Production pH=7.00' !!
				cc = i+1
				while cc < len(arg):
					if arg[cc] not in flags:
						a_t = arg[cc]
					else:
						i = cc - 1
						break

					b_t = ''
					cc_t = cc +1
					while cc_t < len(arg):
						if arg[cc_t] not in flags and '.dat' not in arg[cc_t]:
							if len(b_t) > 0:
								b_t += ' '
							b_t += arg[cc_t]
						else:
							cc = cc_t -1
							break
						cc_t += 1

					if len(b_t) > 0:
						File.append( (a_t, b_t) )
					else:
						i = cc  
						break
					cc += 1
				continue
			elif arg[i].lower() == "-stitle":
				supertitle = arg[i+1]
				continue
			elif arg[i].lower() == "-lblcrd":
				#temp_ar = arg[i+1][1:-1].split(',')
				#tempx, tempy = arg
				labx = float(arg[i+1])
				laby = float(arg[i+2])
				continue
			elif arg[i].lower() == "-fram2time":
				fram2time = True
				framstp = float(arg[i+1])
				continue
			elif arg[i].lower() == "-nanosec":
				nano = True
				#continue
			elif arg[i].lower() == "-dpi":
				dpi = arg[i+1]
				continue
			cut +=1

		else: #cut!= i means that the current arg[i] was used in the previous iteration
			cut = i+1

	#for i in [7.0,8.0,9.0,10.0]:
	#	File.append( ('Gpu-Ultra_cphmd-petase-water/ph%.2f_%s.dat'%(i,anatp),'pH=%d'%i) ) #PETase_CpHMD_water_analysis/

	if not inst_only and not version_only:
		ob4 = Analysis_plot(type=tpe,names=File, analysisType=anatp, suptitle=supertitle, largerYaxis=True, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi)
