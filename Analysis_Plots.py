import numpy as np
import matplotlib.pyplot as plt

class Analysis_plot:
	_Types = ['one','two','2ph_2merge','four','four_2merge','allin_one']
	def __init__(self, type = 'one', names=[('file.dat','analysis_title')], analysisType = 'rmsd', largerYaxis = False, 
	Yenlarger = 2, frameToTime=False, frameStep = 5*10**4, timeStep = 0.004, nanosec = False, suptitle='Titulo geral', 
	labelpx = 35.0, labelpy = 0.50, dpi = 100, label_color = 'darkblue', merge_legend = '', multi_merge_label_loc = (0.5,0.5)):
		'''Parameters:
		
		type: Plot 'one' dataset, 'two'/'2ph_2merge' datasets beside each other, 'four'/'four_2merge' datasets in a matrix fashion or multiple datasets in one XY frame ('allin_one').

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
		
		merge_legend: (type=four_2merge) Inset legend. 

		dpi: Sets the picture quality.'''

		self.X                 = []
		self.Y                 = []
		self.frameToTime       = frameToTime
		self.frameStep         = frameStep
		self.timeStep          = timeStep
		self.nanosec           = nanosec
		self.dpi               = dpi
		self.suptitle          = suptitle
		self.label_color       = label_color
		self.labelpx           = labelpx
		self.labelpy           = labelpy
		self.merge_legend      = merge_legend
		self.restriction_break = False
		if type == 'allin_one':
			self.restriction_break = True

		if 'radgyr' in analysisType:
			self.ana_type = 'Raio de giro (Angstrom)'
		elif 'dist' in analysisType:
			self.ana_type = 'Distância (Angstrom)'
		else:
			self.ana_type = analysisType.upper()+' (Angstrom)'
		XTlabels = self.XY(files=names)
		if XTlabels != -1 and type == 'one':
			self.plot_one(X=self.X[0], Y=self.Y[0], Xaxis=XTlabels[0][0], 
			name= XTlabels[0][1], EnlargeYaxis=largerYaxis, Ymax=Yenlarger)
		elif XTlabels != -1 and type == 'two':
			self.plot_two(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[0][1], XTlabels[1][1]])
		elif XTlabels != -1 and type == '2ph_2merge': 
			self.plot_two_2merge(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[i][1] for i in range(len(XTlabels))])
		elif XTlabels != -1 and type == 'four':
			self.plot_four(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[i][1] for i in range(len(XTlabels))]) #[XTlabels[0][1], XTlabels[1][1], XTlabels[2][1], XTlabels[3][1]]
		elif XTlabels != -1 and type == 'four_2merge': 
			self.plot_four_2merge(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
			name= [XTlabels[i][1] for i in range(len(XTlabels))])
		elif XTlabels != -1 and type == 'allin_one':
			self.multi_in_one(labels=[XTlabels[i][1] for i in range(len(XTlabels))], label_loc=multi_merge_label_loc, titulo =suptitle, eixoy=self.ana_type, eixox=XTlabels[0][0])
		else:
			print("Type must in the list:",self._Types, "!\n")

	def XY(self, files = [('file1.dat','title1'), ('file2.dat','title2')]):
		''' Gets X and Y datasets and returns its respectives x-label and title as a list of tuples.'''

		if len(files) in [1,2,4,8] or self.restriction_break:
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
			axes = plt.axes()
			axes.set_xlim([0,X[len(X)-1]])
			axes.set_ylim([0,Ymax*Y[len(Y)-1]])
		plt.plot(X,Y)
		plt.title(name)
		plt.ylabel(self.ana_type)
		plt.xlabel(Xaxis)
		plt.grid(True)
		plt.show()

	def plot_two(self, X = [[], []], Y = [[], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, dpi=self.dpi)
		ax1.plot(X[0],Y[0])
		ax1.set_title(name[0])
		ax1.grid()
		ax2.plot(X[1],Y[1])
		ax2.set_title(name[1])
		ax2.grid()
		for ts in [ax1,ax2]:
			ts.set(xlabel=Xaxis, ylabel=self.ana_type)
		#pyplot.subplots can hide redundant axes
		for ts in [ax1,ax2]:
			ts.label_outer()
		plt.show()

	def plot_two_2merge(self, X = [[], []], Y = [[], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, dpi=self.dpi)
		big_pp = 0
		axs = [ax1, ax2]
		for ap in axs:
			ap.plot(X[big_pp],Y[big_pp], label= self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= self.merge_legend[1])
			ap.legend()
			ap.set_title(name[big_pp])
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid()
			big_pp += 2
		
		for ts in axs:
			ts.set(xlabel=Xaxis, ylabel=self.ana_type)
		#pyplot.subplots can hide redundant axes
		for ts in axs:
			ts.label_outer()

		plt.suptitle(self.suptitle)
		plt.show()

	def plot_four(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		plt.figure(dpi=self.dpi)
		ax1 = plt.subplot(221)
		plt.plot(X[0],Y[0])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[0], color=self.label_color)
		plt.ylabel(self.ana_type)
		plt.grid(True)
		ax2 = plt.subplot(222, sharey=ax1)		
		plt.plot(X[1],Y[1])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[1], color=self.label_color)
		plt.grid(True)
		ax3 = plt.subplot(223, sharex=ax1)
		plt.plot(X[2],Y[2])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[2], color=self.label_color)
		plt.ylabel(self.ana_type)
		plt.xlabel(Xaxis)
		plt.grid(True)
		ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
		plt.plot(X[3],Y[3])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[3], color=self.label_color)
		plt.xlabel(Xaxis)
		plt.grid(True)
		plt.suptitle(self.suptitle)
		plt.show()

	def plot_four_2merge(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, dpi=self.dpi)
		big_pp = 0
		axs = [ax1, ax2, ax3, ax4]
		for ap in axs:
			ap.plot(X[big_pp],Y[big_pp], label= self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= self.merge_legend[1])
			ap.legend()
			ap.set_title(name[big_pp])
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid()
			big_pp += 2

		for ts in axs:
			ts.set(xlabel=Xaxis, ylabel=self.ana_type)
		#pyplot.subplots can hide redundant axes
		for ts in axs:
			ts.label_outer()

		plt.suptitle(self.suptitle)
		plt.show()

	def multi_in_one(self, labels=['Run 0','Replicata 1','...'], label_loc=(0.5,0.5), titulo ="Ajuste da curva T1", eixoy="Y(X)", eixox="X"):
    	# X = [[],[],'...'], Y = [[],[],'...']
		# X := self.X ;      Y := self.Y
		if len(self.X) != len(self.Y):
			print('For some reason the number of Y datasets is different from de number of X datasets.')
			return -1

		plt.figure(dpi=self.dpi)

		for i in range(len(self.X)):
			plt.plot(self.X[i], self.Y[i])

		plt.legend(labels,loc=label_loc)
		plt.title(titulo)
		plt.ylabel(eixoy)
		plt.xlabel(eixox)
		plt.show()

if __name__ == "__main__":
	import sys
	arg = sys.argv
	version_n = 1.1
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




Copyright (C) 2021-2022  Braga, B. C. 
e-mail: bruno.braga@ufms.br 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

'''
	analysis_name = ['rmsd','rmsf','radgyr','dist']
	dic_Type      = {1:'one', 13:'allin_one', 2:'two', 22: '2ph_2merge', 4:'four', 42: 'four_2merge'}

	#default keys
	version_only = False
	inst_only    = False
	anatp        = analysis_name[3]
	tpe          = dic_Type[13]
	color        = 'black' #'red' #'darkblue' #label only
	mut          = ['HIS 237 - GLU','ASP 206 - GLU'][1]
	supertitle   = '' #['Produção - Mutação %s'%mut,'Petase nativa - Produção'][1]
	bckbne_comp  = False 
	merge_legend = ''
	labx, laby  = (100.5, 0.51) #rmsd
	fram2time   = True  # False
	framstp     = 62500 # Ultra #gpu-high: 50 000
	nano        = True  # False
	dpi         = 100

	# special cases
	rareplot        = False
	multi_label_loc = (0.5,0.5)
	

	if tpe == 'four':
		path = ['ASP206_GLU_Hotfix1.1/', 'nativaHotfix1.1/All_backbone/'][1]
		File = []
		for value in ['7.00','8.00','9.00','10.00']:
			File.append( ('%sph%s_%s.dat'%(path,value,anatp),'pH='+value) )
	elif tpe == 'four_2merge':
		path         = ['nativaHotfix1.1/All_backbone/', 'ASP206_GLU_Hotfix1.1/']
		merge_legend = ['Nativa', 'ASP206GLU']
		File         = []
		supertitle   = 'Produção'
		for value in ['7.00','8.00','9.00','10.00']:
			for p in path:
				File.append( ('%sph%s_%s.dat'%(p,value,anatp),'pH='+value) )
	elif tpe == '2ph_2merge':
		oneTwo       = 1
		path         = [['PETase_Dynamics/gpu-ultra/Jan_no_warp/','PETase_Dynamics/gpu-ultra/Julho_no_warpReplicata/'][oneTwo], 'PETase_Dynamics/gpu-ultra/Mutation_Hotfix_1.1/ASP206_GLU_Hotfix1.1/']
		merge_legend = [['NativaS1','NativaS2'][oneTwo], 'D206E']
		File         = []
		supertitle   = '' #'Produção'
		for value in ['7.00','9.00']:
			for p in path:
				File.append( ('%sph%s_%s.dat'%(p,value,anatp),'pH='+value) )
	elif tpe == 'two':
		if bckbne_comp:
			path    = 'nativaHotfix1.1/'
			ph_comp = '9.00'
			File    = [('%sAll_backbone/ph%s_%s.dat'%(path,ph_comp,anatp),'Todo Backbone a pH=%s'%ph_comp),
			('%sCA_C_N/ph%s_%s.dat'%(path,ph_comp,anatp),'CA,C,N do Backbone a pH=%s'%ph_comp)]
		else:
			path = 'PETase_Dynamics/gpu-ultra/Julho_no_warpReplicata/' #'HIS237_GLU/' #'ASP206_GLU-replicata/'
			File = []
			for value in ['7.00','9.00']:
				File.append( ('%sph%s_%s.dat'%(path,value,anatp),'pH='+value) )
	elif tpe == 'allin_one':
		#path            = 'PETase_Dynamics/gpu-ultra/Dock_run/D206E_C8X/MD_only/'
		path2			= 'PETase_Dynamics/gpu-ultra/Dock_run/Nat_C8X/'
		path_list		= ['Run0_MD/','Replicata1_MD/','Replicata2_MD/']
		rareplot        = True
		path2_list		= ['Run0_CpHMD/','Replicata1_CpHMD/','Replicata2_CpHMD/']
		data_file       = 'Ehyd-PETcarb_7.00.dat'
		supertitle      = ''
		multi_label_loc = (0.75,0.75)
		
		File            = []
		File2           = []
		dyn_value       = 1
		for d in path_list: #[:-2]:
			File.append( (path2+d+data_file,'MD %d'%dyn_value) )
			dyn_value += 1
		
		dyn_value       = 1
		for d in path2_list: #[:-1]:
			File2.append( (path2+d+data_file,'CpHMD %d'%dyn_value) )
			dyn_value += 1
		
		'''framstp  = int(framstp/2) # for nativa MD
		for d in ['Zeroth/','Replicata1/','Replicata2/']:
			File.append( (path+d+data_file,'D %d'%dyn_value) )
			dyn_value += 1'''
		
	else:
		path = 'PETase_Dynamics/gpu-ultra/Dock_run/Nat_C8X/'
		rep_type = ['MD','CpHMD'][1]
		#framstp  = int(framstp/2)
		rep  = ['Run0_%s/'%rep_type,'Replicata1_%s/'%rep_type,'Replicata2_%s/'%rep_type]
		File = [('%sEhyd-PETcarb_%s.dat'%(path+rep[2],'7.00'),'')]

	if anatp == 'rmsf':
		fram2time  = False
		nano       = False
		labx, laby = (150, 2.25) 
	elif anatp == 'radgyr':
		labx, laby = (130, 16.8)

	flags = ["&","-v","--version","-h","--help",'-type', '-i', '-anatp','-stitle','-fram2time', '-framstp','-nanosec','-dpi','-lblcrd','-mlbpos']

	cut = 0 # counter for input flags

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
				print("Copyright (C) 2021  Braga, B. C.\nThis program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; use option '-v' for details.\n")
				print("\nUsage:\n\t-v or --version\t\tprints current version, the python libraries needed and the data format expected.\n")
				print("\t-h or --help\t\tprints this message.\n")
				print("\t-type\tQuantity of files used, for data comparison.\n\t\tone: Normal plot with one data file.\n\t\ttwo: Plots two data files with the same axis information (ex:RMSD for pH7 and pH8).\n\t\tfour: Plots four data files with the same axis information (Eg.:RMSD for pH5, pH6, pH7, pH8).\n\t\tallin_one: Plots every dataset given in one XY frame.\n")
				print("\t-anatp\t\tAnalysis type:\n\t\trmsd;\n\t\trmsf;\n\t\tradgyr;\n\t\tdist.\n")
				print("\t-i\t\tinput data file(s) with a name for the plot (separated by space).\n\t\t\tEg.: -i ph7.00_rmsd.dat pH=7.00\n")
				print("\t-stitle\t\t(Valid only for type 'four') Title for comparison plot.\n")
				print("\t-lblcrd\t(Valid only for type 'four') Label coords.\n")
				print("\t-mlbpos\t(Valid only for type 'allin_one') Label position based on axis percentage. Eg.: -mlbpos 0.5 0.8\n")
				print("\t-fram2time\t(For anatp = rmsd or radgyr) Sets X-axis will change from frame to time (pico seconds). Must inform frame-step conversion.\n\t\t\tEg.: -fram2time 10000\n")
				print("\t-nanosec\t\t(For anatp = rmsd or radgyr) Sets time intervals on X-axis to nanoseconds.\n")
				print("\t-dpi\t\tSets the plot quality (recommended to use only with type=one). Default: 100\n")
				print("\nExamples:\n\t$ python3 Analysis_Plots.py -type one -anatp rmsd -i ph7.00_rmsd.dat Production pH=7.00 -fram2time 10000 -dpi 120\n")
				print("\n\t$ python3 Analysis_Plots.py -type two -anatp rmsf -i ph7.00_rmsf.dat Production pH=7.00 ph8.00_rmsf.dat Production pH=8.00\n")
				print("\n\t$ python3 Analysis_Plots.py -type four -stitle Production -anatp radgyr -lblcrd (16.61,200) -i ph7.00_radgyr.dat pH=7.00 ph8.00_radgyr.dat pH=8.00 ph9.00_radgyr.dat pH=9.00 ph10.00_radgyr.dat pH=10.00 -fram2time 10000 -nanosec\n")
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
				sup_t = ''
				cc = i + 1
				while cc < len(arg):
					if arg[cc] not in flags:
						if len(sup_t) > 0:
							sup_t += ' '
						sup_t += arg[cc]
					else:
						i = cc -1
						break
					cc += 1
				supertitle = sup_t
				continue
			elif arg[i].lower() == "-lblcrd":
				#temp_ar = arg[i+1][1:-1].split(',')
				#tempx, tempy = arg
				labx = float(arg[i+1])
				laby = float(arg[i+2])
				continue
			elif arg[i].lower() == "-mlbpos":
				multi_label_loc = (float(arg[i+1]),float(arg[i+2]))
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

	if not inst_only and not version_only:
		if not rareplot:
			ob4 = Analysis_plot(type=tpe,names=File, analysisType=anatp, suptitle=supertitle, largerYaxis=True, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi, label_color=color, merge_legend=merge_legend, multi_merge_label_loc=multi_label_loc )
		else:
			# creating an empty object
			ob4   = Analysis_plot(type='tpe',names=File2, analysisType=anatp, suptitle=supertitle, largerYaxis=True, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi, label_color=color, merge_legend=merge_legend, multi_merge_label_loc=multi_label_loc )
			ob4.X = []
			ob4.Y = []
			
			
			ob4.restriction_break = True
			CpHMDXYlabels         = ob4.XY(files=File2)
			#print(File2) #ok
			print(CpHMDXYlabels)
			labels_cphmd   = [CpHMDXYlabels[i][1] for i in range(len(CpHMDXYlabels))]
			CpHMDeixox     = CpHMDXYlabels[0][0]
			Xcphmd         = ob4.X
			Ycphmd         = ob4.Y
			
			ob4.X          = []
			ob4.Y          = []
			ob4.frameStep = int(framstp/2)
			MDXYlabels = ob4.XY(files=File)
			labels_md  = [MDXYlabels[i][1] for i in range(len(MDXYlabels))]
			MDeixox    = MDXYlabels[0][0]
			
			dual_labels = []
			dual_eixox  = MDeixox # it's the same as CpHMDeixox 
			dual_labels.extend( labels_md )
			dual_labels.extend( labels_cphmd )
			ob4.X.extend( Xcphmd )
			ob4.Y.extend( Ycphmd )
			
			if MDXYlabels != -1 and CpHMDXYlabels != 1:
				ob4.multi_in_one(labels=dual_labels, label_loc=multi_label_loc, titulo =supertitle, eixoy=ob4.ana_type, eixox=dual_eixox)
			else:
				print("Sheesh!")
