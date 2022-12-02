import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import re

class Analysis_plot:
	_Types = ['one','two','2ph_2merge','1c3p','four','four_2merge','allin_one']
	def __init__(self, type = 'one', print_mean=False, names=[('file.dat','analysis_title')], analysisType = 'rmsd',mult_ana_plot=[],
	frameToTime=False, frameStep = 5*10**4, timeStep = 0.004, nanosec = False, suptitle='Titulo geral', 
	labelpx = 35.0, labelpy = 0.50, dpi = 100, label_color = 'darkblue', merge_legend = '', multi_merge_label_loc = (0.5,0.5),
	forced_3Dzaxis=['ph','index'][1], forced_mean=False, fmv=4.75, eng=False, fontsize=10, mmpbsa=False,
	mmpbsa_inset=[False,True][0], mmpbsa_inset_range=(169,210), mmpbsa_inset_tick=4, mmpbsa_inset_XY=(0.1,0.125),
	plot_interface=False, mmpbsa_cut=-0.5, halving_ids=[], debug=False, file_name='Figure.jpeg'):
		'''Parameters:
		
		mult_ana_plot: Default: Empty ([]). If not empty it will have the same size as 'names' and it will correspond to the analysis on each file in 'names'.

		type: Plot 'one' dataset, 'two'/'2ph_2merge' datasets beside each other, 'four'/'four_2merge' datasets in a matrix fashion or multiple datasets in one XY frame ('allin_one').

		print_mean: (For type 'allin_one') prints the mean line of all data

		names: Data Files, the format expected is .dat, if any other is used do not expect a pretty result.
		
		analysisType: Y-axis label between: rmsd; rmsf; radgyr; dist; edcomp; prot_state.
		
		frameToTime: Boolean argument. For rmsd and radgyr analysis. If True X-axis will change from frame to time (pico seconds).
			
		frameStep: Integer value. Amount of 'timeStep's used in the production. Default 50 000 (gpu-high mode from ASM.py).
		
		timeStep: 0.002 ps or 0.004 ps for HMR. Default 0.004 ps.

		nanosec: Boolean argument. If True "time" marks on X-axis will be on nanoseconds.
		
		suptitle: Title of the whole plot. Giving titles to every subplot when type is 'four' decreases the pict. quality.
		
		labelpx: X position of internal labelling.
		
		labelpy: Y position of internal labelling.
		
		merge_legend: (type=four_2merge) Inset legend. 

		dpi: Sets the picture quality.
		
		eng: Boolean argument.
				If 'True' sets default texts on english.
		 		If 'False' sets default texts on portuguese.
		mmpbsa: Boolean argument. If True the data will be considered the 'decode_mmpbsa.py' format.
		mmpbsa_cut: Energy limit for which residue is highlighted on the plot. 
		'''
		self.file_name = file_name
		self.max_state = 0
		self.ID_shift = 28
		self.lang_set          = 0
		if eng:
			self.lang_set      = 1
		# Set the font size. Either an relative value of
		# 'xx-small', 'x-small', 'small', 'medium', 'large', 
		# 'x-large', 'xx-large' or an absolute 
		# font size, e.g., 12.
		self.fontsize          = fontsize
		# weight: A numeric value in the range 0-1000 or one of
		# 'ultralight', 'light', 'normal' (default), 'regular', 
		# 'book', 'medium', 'roman', 'semibold', 'demibold', 'demi',
		# 'bold', 'heavy', 'extra bold', 'black'
		# style: Either 'normal' (default), 'italic' or 'oblique'.
		self.forced_mean       = forced_mean
		self.forced_mean_value = fmv
		self.X                 = []
		self.Y                 = []
		self.halving_ids       = halving_ids
		self.frameToTime       = frameToTime
		self.frameStep         = frameStep
		self.timeStep          = timeStep
		self.nanosec           = nanosec
		self.dpi               = dpi
		self.suptitle          = suptitle
		self.supertitle        = False
		if self.suptitle != '':
			self.supertitle = True
		self.label_color       = label_color
		self.labelpx           = labelpx
		self.labelpy           = labelpy
		self.merge_legend      = merge_legend
		self.restriction_break = False
		self.mmpbsa               = mmpbsa
		self.mmpbsa_inset         = mmpbsa_inset
		self.mmpbsa_inset_XY      = mmpbsa_inset_XY
		self.mmpbsa_inset_restick = mmpbsa_inset_tick
		self.mmpbsa_inset_range   = range(mmpbsa_inset_range[0],mmpbsa_inset_range[1])
		self.mmpbsa_cut           = mmpbsa_cut
		self.mmpbsa_res           = []
		self.plot_interface    = plot_interface

		if type == 'allin_one' or type == 'allin_one_f3d' or type == '1c3p':
			self.restriction_break = True

		if 'radgyr' in analysisType:
			self.ana_type = ['Raio de giro ($\AA$)','Radgyr ($\AA$)'][self.lang_set]
		elif 'dist' in analysisType:
			self.ana_type = ['Distância ($\AA$)','Distance ($\AA$)'][self.lang_set]
		elif 'prot_state' in analysisType:
			self.ana_type = ['Estado de Protonação','Protonation State'][self.lang_set]
		else:
			self.ana_type = analysisType.upper()+' ($\AA$)'
		
		self.mult_ana = False
		self.ana_list = []
		if len(mult_ana_plot) == len(names):
			self.mult_ana = True
			self.ana_list = mult_ana_plot
			##
			if type in ['2ph_2merge','four_2merge','allin_one','allin_one_f3d']:
				XTlabels = -1
				print("Can not merge different types of analysis into one plot.")
			else:
				XTlabels=[]
				new_list = []
				j=-1
				for i in self.ana_list:
					j +=1
					if i =='edcomp':
						new_list.append(['Decomposição de energia\n(kcal/mol)','Energy decomposition\n(kcal/mol)'][self.lang_set])
						XTlabels.extend( self.XY_decomp(files=names) )
					else:
						if i=='radgyr':
							temp = ['Raio de giro ($\AA$)','Radgyr ($\AA$)'][self.lang_set]
						elif i=='dist':
							temp = ['Distância ($\AA$)','Distance ($\AA$)'][self.lang_set] 
						elif i=='prot_state':
							temp = ['Estado de Protonação','Protonation State'][self.lang_set] 
						else:
							temp = analysisType.upper()+' ($\AA$)'
						new_list.append( temp )
						if debug:
							print(temp)
						XTlabels.extend( self.XY(files=[names[j]],ana_type=temp) )
				self.ana_list = new_list
		else:
			if len(mult_ana_plot) >0 and len(names)> 2:
				print("'-multana' doesn't include the analysis for all the files!")
				print("-multana: ",mult_ana_plot)
				print("-i:", names)
				XTlabels = -1
			elif self.mmpbsa:
				XTlabels = self.XY_decomp(files=names)
			else:
				XTlabels = self.XY(files=names,ana_type=self.ana_type)

		if debug:
			print(XTlabels, '\nnames:',names)
			#print(self.X)
			#print(self.Y)
			print("tamanho X e Y:",len(self.X[0]),len(self.Y[0]))
			
		else:
			if XTlabels != -1 and type == 'one':
				self.plot_one(X=self.X[0], Y=self.Y[0], Xaxis=XTlabels[0][0], 
				name= XTlabels[0][1])
			elif XTlabels != -1 and type == 'two':
				self.plot_two(X=self.X, Y=self.Y, Xaxis=[XTlabels[0][0],XTlabels[1][0]],
				name= [XTlabels[0][1], XTlabels[1][1]])
			elif XTlabels != -1 and type == '1c3p':
				self.plot_1c3p(X=self.X, Y=self.Y, Xaxis=[XTlabels[0][0],XTlabels[1][0],XTlabels[2][0]],
				name= [XTlabels[0][1], XTlabels[1][1], XTlabels[2][1]])
			elif XTlabels != -1 and type == 'four':
				self.plot_four(X=self.X, Y=self.Y, Xaxis=[XTlabels[0][0],XTlabels[1][0],XTlabels[2][0],XTlabels[3][0]],
				name= [XTlabels[i][1] for i in range(len(XTlabels))]) #[XTlabels[0][1], XTlabels[1][1], XTlabels[2][1], XTlabels[3][1]]
			elif XTlabels != -1 and type == '2ph_2merge': 
				self.plot_two_2merge(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
				name= [XTlabels[i][1] for i in range(len(XTlabels))])
			elif XTlabels != -1 and type == 'four_2merge': 
				self.plot_four_2merge(X=self.X, Y=self.Y, Xaxis=XTlabels[0][0],
				name= [XTlabels[i][1] for i in range(len(XTlabels))])
			elif XTlabels != -1 and type == 'allin_one':
				self.multi_in_one(labels=[XTlabels[i][1] for i in range(len(XTlabels))], label_loc=multi_merge_label_loc, eixoy=self.ana_type, eixox=XTlabels[0][0],mean=print_mean)
			elif XTlabels != -1 and type == 'allin_one_f3d':
				self.multi_in_one_forced_3d(Zaxis=forced_3Dzaxis,labels=[XTlabels[i][1] for i in range(len(XTlabels))], label_loc=multi_merge_label_loc, titulo =suptitle, eixoy=self.ana_type, eixox=XTlabels[0][0],mean=print_mean)
			else:
				print("Type must in the list:",self._Types, "!\n")

	def XY(self, files = [('file1.dat','title1'), ('file2.dat','title2')],ana_type='self.ana_type'):
		''' Gets X and Y datasets and returns its respectives x-label and title as a list of tuples.'''
		prot_state= False
		if ['Estado de Protonação','Protonation State'][self.lang_set] in ana_type:
			prot_state = True

		if len(files) in [1,2,4,8] or self.restriction_break:
			xname = 'T'
			if prot_state:
				xname = ['Tempo (ns)','Time (ns)'][self.lang_set]
			XTlabel = []
			for fs in files:
				f = open(fs[0])
				t = f.readlines()
				f.close()
				plot_title = fs[1]
				mdid = int(re.findall(r'[0-9]+',plot_title)[-1])
				if mdid in self.halving_ids:
					half_data = True
				else:
					half_data = False
				x = []
				y = []
				xcount = 1 # specific for cutting data
				for i in t:
					data = i.split()
					if i != '' and i[0] == '#':
						if 'RMSF' in ana_type:
							xname = ['Resíduo','Residue'][self.lang_set]
						elif not prot_state and not self.frameToTime:
							xname = data[0][1:]
							if 'frame' in xname.lower():
								xname = ['Instância','Frame'][self.lang_set]
						elif not prot_state:
							if self.nanosec:
								xname = ['Tempo (ns)','Time (ns)'][self.lang_set]
							else:
								xname = ['Tempo (ps)','Time (ps)'][self.lang_set]
					elif len(data) == 2:
						if prot_state:
							x.append( float(data[1]) )
							y.append( int(data[0]) )
							if int(data[0]) > self.max_state:
								self.max_state = int(data[0])
						else:
							x_value = float(data[0])
							y.append( float(data[1]) )
							if half_data:
								if x_value%2 ==0:
									continue
								x_value = xcount

						if 'RMSF' in ana_type:
							x.append( int(x_value) )
						elif self.frameToTime and not prot_state:
							if self.nanosec:
								x.append( x_value* self.frameStep * self.timeStep/1000)
							else:
								x.append( x_value* self.frameStep * self.timeStep)
						elif not prot_state:
							x.append( x_value )
						xcount += 1
				self.X.append(x)
				self.Y.append(y)
				XTlabel.append( (xname,plot_title) )
			return XTlabel
		else:
			return -1

	def XY_decomp(self, files = [('file1.dat','title1'), ('file2.dat','title2')]):
		''' Specific for 'decode_mmpbsa.py' output format.'''

		if len(files) in [1,2,4,8] or self.restriction_break:
			xname = ['Resíduo','Residue'][self.lang_set]
			if not self.mult_ana:
				self.ana_type = ['Decomposição de energia\n(kcal/mol)','Energy decomposition\n(kcal/mol)'][self.lang_set]
			XTlabel = []
			for fs in files:
				f = open(fs[0])
				t = f.readlines()
				f.close()
				plot_title = fs[1]
				x = []
				y = []
				if self.restriction_break:
					decomp_res = []
				default_counter = {'HIP': 'HIS', 'AS4': 'ASP', 'GL4': 'GLU'}
				for i in t:
					data = i.split()
					if "Total" not in i and i != '' and i != '\n':
						if 'R' in data[2]:
							if data[0] in default_counter:
								if self.restriction_break:
									decomp_res.append( default_counter[data[0]] )
								else:
									self.mmpbsa_res.append( default_counter[data[0]] )
							else:
								if self.restriction_break:
									decomp_res.append( data[0] )
								else:
									self.mmpbsa_res.append( data[0] )
							x.append( int(data[1]) )
							y.append( float(data[3]) )
				if self.restriction_break:
					self.mmpbsa_res.append( decomp_res )
				self.X.append(x)
				self.Y.append(y)
				XTlabel.append( (xname,plot_title) )
				#looks unnecessary but I want the same format of the other analysis!
			return XTlabel
		else:
			return -1

	def peaks(self, X=[], Y=[]):
		notes      = []
		false_peak = []
		for i in range(len(Y)):
			if Y[i] <= self.mmpbsa_cut or Y[i] >= -1*self.mmpbsa_cut:
				false_peak.append( (X[i],Y[i]) )
				if len(false_peak)==1:
					notes.append( (X[i],Y[i]) )
				elif false_peak[-1][0]-false_peak[-2][0] > 1:
					#i don't want to see a bunch of neighbouring residues 
					notes.append( (X[i],Y[i]) )
		return notes	

	def plot_one(self, X = [], Y = [], Xaxis = "Frames", name = "Título"):
		'''Plots one X-Y dataset.'''

		mark_color="lightsteelblue"
		plot_color="cornflowerblue"
		inset_color = "darkslategray"
		
		if not self.mmpbsa: 
			fig = plt.figure(dpi=self.dpi)
			#fig, ax1 = plt.subplots(nrows=1, ncols=1, dpi=self.dpi)
			if ['Estado de Protonação','Protonation State'][self.lang_set] in self.ana_type:
				plt.plot(X,Y,'o',ms=1)
				plt.yticks(range(0,self.max_state+1))
			else:
				plt.plot(X,Y)
			plt.title(name, fontsize=self.fontsize)
			plt.ylabel(self.ana_type, fontsize=self.fontsize)
			plt.xlabel(Xaxis, fontsize=self.fontsize)
		else:
			fig, ax1 = plt.subplots(nrows=1, ncols=1, dpi=self.dpi)
			ax1.plot(X,Y,color=plot_color)
			ax1.set_ylabel(self.ana_type, fontsize=self.fontsize)
			ax1.set_xlabel(Xaxis, fontsize=self.fontsize)
			ax1.grid(True)
			ax1.set_ylim(min(Y)-2,max(Y)+0.25)
			# inset data res 170-210
			if self.mmpbsa_inset:
				inset1_x   = []
				inset1_y   = []
				for i in self.mmpbsa_inset_range:
					inset1_x.append( X[i] )
					inset1_y.append( Y[i] )
				inset1 = plt.axes([0,0,1,1],label='upinset')
				inset1.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset1.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip1 = InsetPosition(ax1, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				inset1.set_axes_locator(ip1)
				mark_inset(ax1, inset1, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset1.plot(inset1_x,inset1_y,color=plot_color)
				inset1.axvline(x=178,color='green')
				inset1.axvline(x=209,color='red')
			
			note = self.peaks(X=X,Y=Y)
			for i in note:
				text_color  ='black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.4][2]
				if i[0] == 180:
					cor = [0.2,0.4][1]
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset1.text(i[0],i[1],self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize)
					else:
						ax1.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
				else:		
					ax1.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_two(self, X = [[], []], Y = [[], []], Xaxis = ["Frames"], name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''
		
		mark_color  ="lightsteelblue"
		plot_color  ="cornflowerblue"
		inset_color = "darkslategray"
		share=True
		if self.mult_ana:
			share = False
		fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=share, dpi=self.dpi)
		anatest = ['Estado de Protonação','Protonation State'][self.lang_set]
		if self.mult_ana:
			jj = -1
			for ts in [ax1, ax2]:
				jj+=1 
				if self.ana_list[0]==anatest:
					plt.setp(ts, yticks=range(0,self.max_state+1))
				'''else:
					plt.setp(ts, yticks=y_tick)'''

		if self.mult_ana and self.ana_list[0]==anatest or anatest in self.ana_type:
			ax1.plot(X[0],Y[0],'o',ms=1)
		else:
			ax1.plot(X[0],Y[0],color=plot_color)
			ax1.set_ylim(min(Y[0])-0.5,max(Y[0])+0.5)
		ax1.set_title(name[0], fontsize=self.fontsize)
		ax1.grid()

		if self.mult_ana and self.ana_list[0]=='edcomp' or self.mmpbsa:
			# inset data res 170-210
			ax1.set_xlabel(Xaxis[0],fontsize=self.fontsize)
			if self.mult_ana:
				ylb = self.ana_list[0]
			else:
				ylb = self.ana_type
			ax1.set_ylabel(ylb,fontsize=self.fontsize)
			if self.mmpbsa_inset:
				inset1_x   = []
				inset1_y   = []
				for i in self.mmpbsa_inset_range:
					inset1_x.append( X[0][i] )
					inset1_y.append( Y[0][i] )
				inset1 = plt.axes([0,0,1,1],label='upinset')
				inset1.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset1.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip1 = InsetPosition(ax1, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				inset1.set_axes_locator(ip1)
				mark_inset(ax1, inset1, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset1.plot(inset1_x,inset1_y,color=plot_color)
				inset1.axvline(x=178,color='green')
				inset1.axvline(x=209,color='red')
			
			note = self.peaks(X=X[0],Y=Y[0])
			for i in note:
				text_color  ='black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.4][2]
				if i[0] == 180:
					cor = [0.2,0.4][1]
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset1.text(i[0],i[1],self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize)
					else:
						ax1.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
				else:		
					ax1.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
		else:
			if self.mult_ana:
				ylb = self.ana_list[0]
			else:
				ylb = self.ana_type
			ax1.set_xlabel(Xaxis[0], fontsize=self.fontsize)
			ax1.set_ylabel(ylb, fontsize=self.fontsize)	
		
		if self.mult_ana and self.ana_list[1]==anatest or anatest in self.ana_type:
			ax2.plot(X[1],Y[1],'o',ms=1)
		else:
			ax2.plot(X[1],Y[1],color=plot_color)
			ax2.set_ylim(min(Y[1])-0.5,max(Y[1])+2)
		ax2.set_title(name[1], fontsize=self.fontsize)
		ax2.grid()
		#ax2.set_ylim([-4.5,0.2])
		if self.mult_ana and self.ana_list[1]=='edcomp' or self.mmpbsa:
			ax2.set_xlabel(Xaxis[1],fontsize=self.fontsize)
			if self.mult_ana:
				ylb = self.ana_list[1]
			else:
				ylb = self.ana_type
			ax2.set_ylabel(ylb,fontsize=self.fontsize)
			
			if self.mmpbsa_inset:
				inset2_x   = []
				inset2_y   = []
				for i in self.mmpbsa_inset_range:
					inset2_x.append( X[1][i] )
					inset2_y.append( Y[1][i] )
				inset2 = plt.axes([0,0,1,1],label='downinset')
				inset2.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset2.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip2 = InsetPosition(ax2, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				mark_inset(ax2, inset2, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset2.set_axes_locator(ip2)
				inset2.plot(inset2_x,inset2_y,color=plot_color)
				inset2.axvline(x=209,color='red')

			note = self.peaks(X=X[1],Y=Y[1])
			for i in note:
				#ajustments with 'cor' are manual for now 
				text_color = 'black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.6][2]
				if i[0] == 180:
					cor = 0.2
				elif i[0] == 59:
					cor = [0,-0.5][1]	
				elif i[0] == 133:
					cor = [0,0.15,1.2][2]	
				elif i[0] == 157:
					cor = [0,0.35][1]	
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset2.text(i[0],i[1]+cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize)
					else:
						ax2.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
				else:
					ax2.text(i[0], i[1]-cor, self.mmpbsa_res[i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize)
		else:
			if self.mult_ana:
				ylb = self.ana_list[1]
			else:
				ylb = self.ana_type
			ax2.set_xlabel(Xaxis[1], fontsize=self.fontsize)
			ax2.set_ylabel(ylb, fontsize=self.fontsize)
			
		#pyplot.subplots can hide redundant axes
		for ts in [ax1,ax2]:
			ts.label_outer()
		
		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_1c3p(self, X = [[], []], Y = [[], []], Xaxis = ["Frames"], name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''
		
		#testar com multi ana!!!!!!!!!!!!!11

		mark_color  ="lightsteelblue"
		plot_color  ="cornflowerblue"
		inset_color = "darkslategray"
		vline_c     = 'purple'
		vline_thickness = [0.25,0.75][1]
		fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True, dpi=self.dpi)
		max_y = -50
		min_y = +50
		for yy in Y:
			max_y = max(max_y,max(yy))
			min_y = min(min_y,min(yy))
			#print(max_y,min_y)
		y_tick=np.arange(round(min_y,1),round(max_y,1)+1,round((max_y-min_y)/4,1))
		
		anatest = ['Estado de Protonação','Protonation State'][self.lang_set]
		anatest2 = ['Decomposição de energia\n(kcal/mol)','Energy decomposition\n(kcal/mol)'][self.lang_set] 
		if not self.mult_ana:
			plt.setp((ax1, ax2, ax3), yticks=y_tick)
		else:
			jj = -1
			for ts in [ax1, ax2, ax3]:
				jj+=1 
				if self.ana_list[jj]==anatest:
					plt.setp(ts, yticks=range(0,self.max_state+1))
				else:
					plt.setp(ts, yticks=y_tick)
		plt.subplots_adjust(hspace=0.3)

		if self.mult_ana and self.ana_list[0]==anatest or anatest in self.ana_type:
			ax1.plot(X[0],Y[0],'o',ms=1)
		elif anatest2 in self.ana_type:
			print('residue correction', self.ID_shift)
			X_0 = []
			for xx in X[0]:
				X_0.append(xx+self.ID_shift)
			ax1.plot(X_0,Y[0],color=plot_color)
			ax1.set_ylim(round(min_y,1)-0.5,round(max_y,1)+1)
		else:
			ax1.plot(X[0],Y[0],color=plot_color)
		ax1.set_title(name[0], fontsize=self.fontsize)
		ax1.grid()
		if self.mult_ana and self.ana_list[0]=='edcomp' or self.mmpbsa:
			ax1.set_xlabel(Xaxis[0],fontsize=self.fontsize)
			if self.mult_ana:
				ylb = self.ana_list[0]
			else:
				ylb = self.ana_type
			ax1.set_ylabel(ylb,fontsize=self.fontsize)
			# inset data res 170-210
			if self.mmpbsa_inset:
				inset1_x   = []
				inset1_y   = []
				for i in self.mmpbsa_inset_range:
					inset1_x.append( X_0[i] )
					inset1_y.append( Y[0][i] )
				inset1 = plt.axes([0,0,1,1],label='upinset')
				inset1.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset1.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip1 = InsetPosition(ax1, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				inset1.set_axes_locator(ip1)
				mark_inset(ax1, inset1, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset1.plot(inset1_x,inset1_y,color=plot_color)
				#inset1.axvline(x=178,color='green')
				#inset1.axvline(x=209,color='red')
			
			note = self.peaks(X=X[0],Y=Y[0])
			for i in note:
				text_color  ='black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.4][2]
				if i[0] == 180:
					cor = [0.2,0.4][1]
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset1.text(i[0]+self.ID_shift,i[1],self.mmpbsa_res[0][i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize*0.5)
					else:
						ax1.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[0][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)
				else:		
					ax1.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[0][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)

		if self.mult_ana and self.ana_list[1]==anatest or anatest in self.ana_type:
			ax2.plot(X[1],Y[1],'o',ms=1)
		elif anatest2 in self.ana_type:
			X_1 = []
			for xx in X[1]:
				X_1.append(xx+self.ID_shift)
			ax2.plot(X_1,Y[1],color=plot_color)
			ax2.set_ylim(round(min_y,1)-0.5,round(max_y,1)+1)
		else:
			ax2.plot(X[1],Y[1],color=plot_color)
		ax2.set_title(name[1], fontsize=self.fontsize)
		ax2.grid()
		#ax2.set_ylim([-4.5,0.2])
		if self.mult_ana and self.ana_list[1]=='edcomp' or self.mmpbsa:
			ax2.set_xlabel(Xaxis[1],fontsize=self.fontsize)
			if self.mult_ana:
				ylb = self.ana_list[1]
			else:
				ylb = self.ana_type
			ax2.set_ylabel(ylb,fontsize=self.fontsize)
			if self.mmpbsa_inset:
				inset2_x   = []
				inset2_y   = []
				for i in self.mmpbsa_inset_range:
					inset2_x.append( X_0[i] )
					inset2_y.append( Y[1][i] )
				inset2 = plt.axes([0,0,1,1],label='downinset')
				inset2.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset2.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip2 = InsetPosition(ax2, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				mark_inset(ax2, inset2, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset2.set_axes_locator(ip2)
				inset2.plot(inset2_x,inset2_y,color=plot_color)
				inset2.axvline(x=209,color='red')

			note = self.peaks(X=X[1],Y=Y[1])
			for i in note:
				#ajustments with 'cor' are manual for now 
				text_color = 'black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.6][0]
				if i[0] == 178:
					text_color = 'green'
					cor = 0.5
				elif i[0] == 132:
					cor = [0,0.15,1.2][2]	
				elif i[0] == 157:
					cor = [0,0.35][1]	
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset2.text(i[0]+self.ID_shift,i[1]+cor, self.mmpbsa_res[1][i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize*0.5)
					else:
						ax2.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[1][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)
				else:
					ax2.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[1][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)

		if self.mult_ana and self.ana_list[2]==anatest or anatest in self.ana_type:
			ax3.plot(X[2],Y[2],'o',ms=1)
		elif anatest2 in self.ana_type:
			X_2 = []
			for xx in X[2]:
				X_2.append(xx+self.ID_shift)
			ax3.plot(X_2,Y[2],color=plot_color)
			ax3.set_ylim(round(min_y,1)-0.5,round(max_y,1)+1)
		else:
			ax3.plot(X[2],Y[2],color=plot_color)
		ax3.set_title(name[2], fontsize=self.fontsize)
		ax3.grid()
		if self.mult_ana and self.ana_list[2]=='edcomp' or self.mmpbsa:
			ax3.set_xlabel(Xaxis[2],fontsize=self.fontsize)
			if self.mult_ana:
				ylb = self.ana_list[2]
			else:
				ylb = self.ana_type
			ax3.set_ylabel(ylb,fontsize=self.fontsize)
			if self.mmpbsa_inset:
				inset3_x   = []
				inset3_y   = []
				for i in self.mmpbsa_inset_range:
					inset3_x.append( X_0[i] )
					inset3_y.append( Y[2][i] )
				inset3 = plt.axes([0,0,1,1],label='downinset')
				inset3.tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
				inset3.set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
				ip3 = InsetPosition(ax3, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
				mark_inset(ax3, inset3, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
				inset3.set_axes_locator(ip3)
				inset3.plot(inset3_x,inset3_y,color=plot_color)
				inset3.axvline(x=209,color='red')

			note = self.peaks(X=X[2],Y=Y[2])
			for i in note:
				#ajustments with 'cor' are manual for now 
				text_color = 'black'
				cor   = 0
				if i[1] >0:
					text_color = 'darkred'
					cor   = [0,0.2,0.6][2]
				if i[0] == 178:
					text_color = 'green'
					cor = 0.8
				elif i[0] == 157:
					cor = [0,0.35][1]	
				if self.mmpbsa_inset:
					if i[0] in self.mmpbsa_inset_range:
						inset3.text(i[0]+self.ID_shift,i[1]+cor, self.mmpbsa_res[2][i[0]-1]+str(i[0]+self.ID_shift), color=inset_color, fontsize=self.fontsize*0.5)
					else:
						ax3.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[2][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)
				else:
					ax3.text(i[0]+self.ID_shift, i[1]-cor, self.mmpbsa_res[2][i[0]-1]+str(i[0]+self.ID_shift), color=text_color, fontsize=self.fontsize*0.5)

		if self.mult_ana:
			counter = -1
			for ts in [ax1,ax2,ax3]:
				counter += 1
				ts.set_xlabel(Xaxis[counter],fontsize=self.fontsize)
				ts.set_ylabel(self.ana_list[counter],fontsize=self.fontsize*0.8)	
		elif self.mmpbsa:
			counter=-1			
			for ts in [ax1,ax2,ax3]:
				counter+=1
				ts.set_xlabel(Xaxis[counter],fontsize=self.fontsize)
				ts.axvline(x=132+self.ID_shift,color=vline_c,lw=vline_thickness,zorder=-1)
				ts.axvline(x=178+self.ID_shift,color=vline_c,lw=vline_thickness,zorder=-1)
				ts.axvline(x=209+self.ID_shift,color=vline_c,lw=vline_thickness,zorder=-1)
			#labelling ax2 only but with bigger font
			ax2.set_ylabel(self.ana_type,fontsize=self.fontsize*1.25)
		else:
			counter = -1
			for ts in [ax1,ax2,ax3]:
				counter += 1
				ts.set_xlabel(Xaxis[counter], fontsize=self.fontsize)
				ts.set_ylabel(self.ana_type, fontsize=self.fontsize)	
		#pyplot.subplots can hide redundant axes
		for ts in [ax1,ax2,ax3]:
			ts.label_outer()
		
		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_two_2merge(self, X = [[], []], Y = [[], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, dpi=self.dpi)
		big_pp = 0
		axs = [ax1, ax2]
		for ap in axs:
			ap.plot(X[big_pp],Y[big_pp], label= self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= self.merge_legend[1])
			ap.set_xlabel(Xaxis,fontsize=self.fontsize)
			ap.set_ylabel(self.ana_type,fontsize=self.fontsize)
			ap.legend(fontsize=self.fontsize)
			ap.set_title(name[big_pp], fontsize=self.fontsize)
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid()
			big_pp += 2
		
		for ts in axs:
			ts.label_outer()

		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_four(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = ["Frames"], name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		plt.figure(dpi=self.dpi)
		ax1 = plt.subplot(221)
		plt.plot(X[0],Y[0])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[0], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[0], fontsize=self.fontsize)
		plt.grid(True)
		ax2 = plt.subplot(222, sharey=ax1)		
		plt.plot(X[1],Y[1])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[1], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[1], fontsize=self.fontsize)
		plt.grid(True)
		ax3 = plt.subplot(223, sharex=ax1)
		plt.plot(X[2],Y[2])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[2], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[2], fontsize=self.fontsize)
		plt.grid(True)
		ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
		plt.plot(X[3],Y[3])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[3], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[3], fontsize=self.fontsize)
		plt.grid(True)
		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		for ts in [ax1,ax2,ax3,ax4]:
			ts.label_outer()

		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_four_2merge(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, dpi=self.dpi)
		big_pp = 0
		axs = [ax1, ax2, ax3, ax4]
		for ap in axs:
			ap.plot(X[big_pp],Y[big_pp], label= self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= self.merge_legend[1])
			ap.legend(fontsize=self.fontsize)
			ap.set_title(name[big_pp], fontsize=self.fontsize)
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid()
			big_pp += 2

		for ts in axs:
			ts.set(xlabel=Xaxis, ylabel=self.ana_type)
		#pyplot.subplots can hide redundant axes
		for ts in axs:
			ts.label_outer()

		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def multi_in_one(self, labels=['Run 0','Replicata 1','...'], label_loc=(0.5,0.5), eixoy="Y(X)", eixox="X",mean=False):
    	# X = [[],[],'...'], Y = [[],[],'...']
		# X := self.X ;      Y := self.Y
		if len(self.X) != len(self.Y):
			print('For some reason the number of Y datasets is different from de number of X datasets.')
			return -1

		plt.figure(dpi=self.dpi)
		plt.grid(True)
		mean_value = []
		for i in range(len(self.X)):
			mean_value.append( sum(self.Y[i])/len(self.Y[i]) )
			plt.plot(self.X[i], self.Y[i])

		if mean:
			mean_multi = sum(mean_value)/len(mean_value)
			mean_line  = [mean_multi]*len(self.X[0])
			plt.plot(self.X[0], mean_line)
			labels.append('M='+str( round(mean_multi,2) ))
		if self.forced_mean:
			plt.plot(self.X[i], [self.forced_mean_value]*len(self.X[i]))
			labels.append('M='+str( round(self.forced_mean_value,2) ))

		plt.legend(labels,loc=label_loc, fontsize=self.fontsize)
		if self.supertitle:
			plt.title(self.suptitle, fontsize=self.fontsize)
		plt.ylabel(eixoy, fontsize=self.fontsize)
		plt.xlabel(eixox, fontsize=self.fontsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def multi_in_one_forced_3d(self,text=['',"(I)","(II)"][0], Zaxis=['ph','index'][1], labels=['Run 0','Replicata 1','...'], label_loc=(0.5,0.5), titulo ="Ajuste da curva T1", eixoy="Y(X)", eixox="X",mean=False):
		if len(self.X) != len(self.Y):
			print('For some reason the number of Y datasets is different from de number of X datasets.')
			return -1

		fig        = plt.figure(dpi=self.dpi)
		ax         = fig.gca(projection='3d')
		mean_value = []
		eixoz      = 'pH'
		if Zaxis.lower() == 'ph':
			change_cc = -1
			temp_ph   = -1
			for i in range(len(self.X)):
				mean_value.append( sum(self.Y[i])/len(self.Y[i]) )
				#labels :: 'pH %c MD %d %c'%(ph,dynNumb,subdivision)
				label_temp  = labels[i].split()
				temp_ph_new = float(label_temp[1])
				if temp_ph != temp_ph_new:
					change_cc += 1
					if mean and change_cc > 0:
						new_ph_temp = mean_value.pop()
						mean_multi  = sum(mean_value)/len(mean_value)
						mean_line   = [mean_multi]*len(self.X[0])
						ax.plot(self.X[0], mean_line, [temp_ph]*len(self.X[i]), label='pH '+str(temp_ph)+' M='+str( round(mean_multi,2) ))
						mean_value  = [new_ph_temp]
					temp_ph = temp_ph_new
				else:
					ph_change = False
				new_label = labels[i][len(label_temp[0]+label_temp[1])+2:]
				ax.plot(self.X[i], self.Y[i], [temp_ph]*len(self.X[i]), label=new_label)
			if mean:
				# for the last pH
				mean_multi = sum(mean_value)/len(mean_value)
				mean_line  = [mean_multi]*len(self.X[0])
				ax.plot(self.X[0], mean_line, [temp_ph]*len(self.X[i]), label='pH '+str(temp_ph)+' M='+str( round(mean_multi,2) ))
		elif Zaxis.lower() == 'index':
			index = 0
			eixoz = eixoy
			eixoy = ['Índice','Index'][self.lang_set]
			for i in range(len(self.X)):
				index += 1
				mean_value.append( sum(self.Y[i])/len(self.Y[i]) )
				#labels :: 'pH %c MD %d %c'%(ph,dynNumb,subdivision)
				ax.plot(self.X[i], [index]*len(self.X[i]), self.Y[i], label=labels[i])
			if mean:
				mean_multi = sum(mean_value)/len(mean_value)
				mean_line  = [mean_multi]*len(self.X[0])
				#index 0
				ax.plot(self.X[0], [1]*len(self.X[i]), mean_line, label='M='+str( round(mean_multi,2) ))
				#ax.plot_surface(self.X[0], [0]*len(self.X[i]), mean_line, label='M='+str(round(mean_multi,2)) )
			if self.forced_mean:
				ax.plot(self.X[0], [1]*len(self.X[i]), [self.forced_mean_value]*len(self.X[0]), label='M='+str( round(self.forced_mean_value,2) ))

		'''legend loc :
Best 0
Upper right 1
Upper left 2
Lower left 3
Lower right 4
Right 5
Center left 6
Center right 7
Lower center 8
Upper center 9
center 10'''
		#frameon makes a box for the legend
		#ax.legend(loc=0,frameon=False,mode="expand")
		ax.legend(bbox_to_anchor=(0, 0.75), borderaxespad=0)

		if self.supertitle:
			fig.suptitle(titulo, fontsize=self.fontsize)
		ax.set_xlabel(eixox, fontsize=self.fontsize)
		ax.set_ylabel(eixoy, fontsize=self.fontsize)
		ax.set_zlabel(eixoz, fontsize=self.fontsize)
		
		if text != '':
			fig.text(0.07, 0.72, text, fontsize=self.fontsize)
		#ax1.plot(Xw,Yw,'o',color=w_color,ms=pointsize)
		
		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

def Default_modifier(on=False,tpe='allin_one', mmpbsa=[True,False][1],
mmpbsa_inset=False, mmpbsa_cut=-0.5, mmpbsa_inset_range=(169,210), interface=True,
anatp='dist',File=[],File2=[],enzyme=['Nat','D206E','D206EH237K'][2],met=[0,1][0],
merge_legend=['NativaS1','NativaS2'],supertitle='Produção',
mean_flag=False,rareplot=False,multi_label_loc=(.74,0.72),
forced_3Dz=['ph','index'][1],forced_mean=False, fmv=4.75,
eng=[True,False][0],fontsize=14,mdlistrange=range(1,3), ph_reverse=False,
igph7md=[1,3,5],igph9md=[1,3,5],igph7cphmd=[],igph9cphmd=[],
res=('LYS209','LYS 237')):

	if not on:
		return -1
	File      = []
	File2     = []
			
	if tpe == 'four' and anatp !='prot_state':
		path = ['ASP206_GLU_Hotfix1.1/', 'nativaHotfix1.1/All_backbone/'][1]
		for value in ['7.00','8.00','9.00','10.00']:
			File.append( ('%sph%s_%s.dat'%(path,value,anatp),'pH='+value) )
	elif tpe == 'four_2merge' and anatp !='prot_state':
		path         = ['nativaHotfix1.1/All_backbone/', 'ASP206_GLU_Hotfix1.1/']
		merge_legend = ['Nativa', 'ASP206GLU']
		for value in ['7.00','8.00','9.00','10.00']:
			for p in path:
				File.append( ('%sph%s_%s.dat'%(p,value,anatp),'pH='+value) )
	elif tpe == '2ph_2merge' and anatp !='prot_state':
		path = 'PETase_Dynamics/gpu-ultra/' 
		s1   = 'Jan_no_warp/'
		s2   = 'Julho_no_warpReplicata/'
		smut = 'Mutation_Hotfix_1.1/ASP206_GLU_Hotfix1.1/'
		for value in ['7.00','9.00']:
			File.append( ('%sph%s_%s.dat'%(path+s1,value,anatp),'pH='+value) )
			File.append( ('%sph%s_%s.dat'%(path+[s2,smut][0],value,anatp),'pH='+value) )
	elif tpe == 'two' and anatp !='prot_state':
		if bckbne_comp:
			path    = 'nativaHotfix1.1/'
			ph_comp = '9.00'
			File    = [('%sAll_backbone/ph%s_%s.dat'%(path,ph_comp,anatp),'Todo Backbone a pH=%s'%ph_comp),
			('%sCA_C_N/ph%s_%s.dat'%(path,ph_comp,anatp),'CA,C,N do Backbone a pH=%s'%ph_comp)]
		else:
			path = 'PETase_Dynamics/gpu-ultra/' 
			s1   = 'Jan_no_warp/'
			s2   = 'Julho_no_warpReplicata/' #'HIS237_GLU/' #'ASP206_GLU-replicata/'
			for value in ['7.00','9.00']:
				File.append( ('%sph%s_%s.dat'%(path+s2,value,anatp),'pH='+value) )
	elif tpe == 'allin_one' or tpe == 'allin_one_f3d' and anatp !='prot_state':
		variant         = {'Nat':0,'D206E':1,'D206EH237K':2,'Nat_nored':4, 'H237K':3}
		path0			= 'PETase_Dynamics/gpu-ultra/Dock_run/'
		path            = [path0+ss for ss in ['Nat_C8X/','D206E_C8X/MD_only/','D206EH237K_C8X/MD_only/','H237K','Nat_nonreduced-C8X/'] ]
		path_list		= [['Run0_MD/','Zeroth/','Run0/','MD1/','MD1/'][variant[enzyme]]]
		path_list.extend( [['Replicata%d_MD/','Replicata%d/','Rep%d/','MD%d','MD%d'][variant[enzyme]]%i for i in mdlistrange] )
		path2_list		= ['Run0_CpHMD/']
		path2_list.extend(['Replicata%d_CpHMD/'%i for i in range(1,3)] )
		data_file       = {'rmsd':['MD_7.00_system.rms'],'rmsf':['MD_7.00_system_rmsf.agr'],'dist':['Ehyd-PETcarb_7.00.dat','Ehyd-PETcarb_9.00.dat']}[anatp]
		if ph_reverse:
			data_file.reverse()
		#linha abaixo para plotar 5 mds da nativa	
		ingore_ph        = ['', ['Ehyd-PETcarb_7.00.dat','Ehyd-PETcarb_9.00.dat'][1]][0]
		ignore_dyn_MD    = {'MD_7.00_system.rms':[],'MD_7.00_system_rmsf.agr':[],'Ehyd-PETcarb_7.00.dat': igph7md,'Ehyd-PETcarb_9.00.dat':igph9md}
		ignore_dyn_CpHMD = {'Ehyd-PETcarb_7.00.dat': igph7cphmd,'Ehyd-PETcarb_9.00.dat':igph9cphmd}
		subdivision_leg  = {'Ehyd-PETcarb_7.00.dat': 'A','Ehyd-PETcarb_9.00.dat':'B'}
		for ph_f in data_file:
			dyn_value    = 1
			for d in path_list:
				if dyn_value not in ignore_dyn_MD[ph_f] and ph_f != ingore_ph:
					# label needed for 'allin_one_f3d' :: 'pH %c MD %d'%(ph_f[-8],dyn_value) 
					if variant[enzyme] == 0 and tpe == 'allin_one' or variant[enzyme] == 4:
						legend = 'MD %d'%(dyn_value)
					else:
						legend = 'pH %c MD %d %c'%(ph_f[-8],dyn_value,subdivision_leg[ph_f])
					File.append( (path[variant[enzyme]]+d+ph_f, legend ) )
				dyn_value += 1
		#File=File2
		if rareplot:
			for ph_f in data_file:
				dyn_value = 1
				for d in path2_list:
					if dyn_value not in ignore_dyn_CpHMD[ph_f]:
						File2.append( (path[3]+d+ph_f,'pH %c CpHMD %d'%(ph_f[-8],dyn_value)) )
					dyn_value += 1
	elif anatp !='prot_state':
		path = ['PETase_Dynamics/gpu-ultra/Dock_run/Nat_C8X/','PETase_Dynamics/gpu-ultra/Dock_run/H237K_C8X/'][0]
		rep_type = ['MD','CpHMD'][0]
		#framstp  = int(framstp/2)
		rep  = ['Run0_%s/'%rep_type,'Replicata1_%s/'%rep_type,'Replicata2_%s/'%rep_type,'MD1/']
		File = [('%sEhyd-PETcarb_%s.dat'%(path+rep[-1],'7.00'),'')]
	
	if mmpbsa and anatp !='prot_state':
		path       = 'IsPETase-BHET_binding/'
		ph         = [7,9][0]
		
		#md_id = (0,0)
		ENZ2id = {'WT-D-DH':(1,4,6),'Nat':(1,5), 'D206E':(4,6), 'D206EH237K':(1,6), 'H237K':(1,1), 'H237Knored':(1,1)}
		md_id = ENZ2id[enzyme]
		metodo = ['MD','CpHMD'][met] 
		#mmpbsa_cut = -1
		if enzyme == 'WT-D-DH':
			File = [('%spH%d-Nat-C8X-MD%d_bind.dat'%(path,ph,md_id[0]),
			'WT-BHET'),
			('%spH%d-D206E-C8X-MD%d_bind.dat'%(path,ph,md_id[1]),
			'D206E-BHET'),
			('%spH%d-D206EH237K-C8X-MD%d_bind.dat'%(path,ph,md_id[2]),
			'D206E/H237K-BHET')]
		elif tpe == 'one':
			File = [('%spH%d-%s-C8X-%s%d_bind.dat'%(path,ph,enzyme,metodo,md_id[0]),
			'%s-BHET pH%d MD%d'%(enzyme,ph,md_id[0]))]
		elif tpe == 'two':
			#coloca uma linha vertical verde no E206 e um roxa no K209!! 
			File = [('%spH%d-%s-C8X-MD%d_bind.dat'%(path,ph,enzyme,md_id[0]),
			'MD%d'%(md_id[0])),
			('%spH%d-%s-C8X-MD%d_bind.dat'%(path,ph,enzyme,md_id[1]),
			'MD%d'%(md_id[1]))]
	elif anatp=='prot_state':
		path       = 'Dev_cpoutAnalyser/'
		ENZ2id = {'Nat':(1,3), 'D206E':(1,6), 'D206EH237K':(1,6)}
		md_id = ENZ2id[enzyme]
		if tpe == 'one':
			File = [('%s_CpH7MD%d-%s_SvT.dat'%(path+enzyme,md_id[1],res[0]),
			'%s'%(res[1]))]
		

	return ( tpe,anatp,File,File2, merge_legend,supertitle, 
	mean_flag,rareplot,multi_label_loc,forced_3Dz,eng,fontsize,
	forced_mean, fmv, mmpbsa, mmpbsa_inset, mmpbsa_inset_range, mmpbsa_cut, interface)

if __name__ == "__main__":
	import sys
	arg = sys.argv
	version_n = 1.2
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
	dic_Type      = {1:'one', 13:'allin_one', 133:'allin_one_f3d', -13:'1c3p', 2:'two', 22: '2ph_2merge', 4:'four', 42: 'four_2merge'}

	#default keys
	interface    = True
	mmpbsa             = False
	mmpbsa_inset       = False
	mmpbsa_inset_range = (169,210)
	mmpbsa_cut         = -0.5
	version_only = False
	inst_only    = False
	anatp        = ['rmsd','rmsf','radgyr','dist','prot_state'][3]
	tpe          = dic_Type[133]
	fontsize     = 10
	eng          = False
	forced_3Dz   = 'index'
	color        = 'black' #'red' #'darkblue' #label only
	mut          = ['HIS 237 - GLU','ASP 206 - GLU'][1]
	supertitle   = '' #['Produção - Mutação %s'%mut,'Petase nativa - Produção'][1]
	bckbne_comp  = False 
	merge_legend = ''
	labx, laby   = (100.5, 0.51) #rmsd
	fram2time    = True  # False
	framstp      = 62500 # Ultra #gpu-high: 50 000
	nano         = True  # False
	dpi          = [100,200,300][2]
	mean_flag    = False
	
	debug        = False # for the arguments only
	debug_inside = [True,False][1] # to help fix plot problems
	debug_dic    = {}

	# special cases #IGNORE THIS
	File         = []
	File2        = []
	merge_legend = ['NativaS1',['NativaS2','D206E'][0]]
	forced_mean  = False
	fmv          = 4.75
	rareplot     = False
	multi_label_loc=(.685,0.5)
	#mmpbsa_cut must be negative!!
	mult_ana_plot = [] # multiple analysis
	file_name = 'Figure.jpeg'

	default_mod = True
	inset_tick      = 10
	mmpbsa_inset_XY = (0.10,0.10)
	halving_ids = [[],[1,2,3]][0]
	met= 0 # MD ## met=1 := CpHMD
	# If the flag '-i' is read then 'default_mod' becames False  
	#tpe=-13 eh 1 coluna 3 linhas
	temp_mod = Default_modifier(on=default_mod,tpe=dic_Type[1],
	anatp=['rmsd','rmsf','radgyr','dist','prot_state'][4],File=File,File2=File2, 
	enzyme=['Nat','D206E','D206EH237K','H237K','H237Knored','Nat_nored','WT-D-DH'][0], mmpbsa=[True,False][1],
	met=1,mmpbsa_cut=[-0.5,-1.5][0], mmpbsa_inset=[True,False][1], mmpbsa_inset_range=[(169,210),(150,160)][0],
	merge_legend=merge_legend,supertitle=['%s-BHET pH7'%['IsPETase','D206EH237K'][1],''][1],
	mean_flag=[True,False][1],rareplot=[True,False][1],interface=[True,False][1],
	multi_label_loc=(.685,0.5),forced_3Dz=['ph','index'][0],
	forced_mean=[True,False][1], fmv=[4.65,4.75,4.97][0], ph_reverse=[True,False][1],
	eng=[True,False][0],fontsize=[8,12,14][1],mdlistrange=range(1,4),
	igph7md=[],igph9md=[1,2,3,4,5,6],igph7cphmd=[],igph9cphmd=[],
	res=[('GL4178','GLU 206'),('LYS209','LYS 237'),('AS4178','ASP 206'),('HIP209','HIS 237'),('TYR59','TYR 87')][3])

	#print("\nmodifiers:",temp_mod)
	
	if anatp == 'rmsf':
		fram2time  = False
		nano       = False
		labx, laby = (150, 2.25) 
	elif anatp == 'radgyr':
		labx, laby = (130, 16.8)

	flags = ["&","-v","--version","-fontsize","-jpeg","-edinsetXYpos","-edinsetX","-edXtick","-edcut","-eng","-h","--help","-type","-index3d","-mean", "-i", "-debug","-anatp","-multana","-stitle","-fram2time", "-framstp","-nanosec","-dpi","-lblcrd","-mlbpos"]

	# Flag verification
	for i in arg:
		if i[0] == '-':
			try:
				if type(float(i)) == type(2.3):
					#Some flags may require a number as the next arg
					continue
			except ValueError:
				if i not in flags:
					print("Unkown Flag used: ", i)
					print("Please check if the CASE is correct as below!")
					print("Flags: ",flags)
					inst_only = True
					break

	if "-jpeg" in arg:
		interface = False

	cut = 0 # counter for input flags
	for i in range(len(arg)):
		if cut == i:
			if arg[i] == "&":
				break
			elif arg[i] == "-v" or arg[i] == "--version":
				version_only = True
				print("Current version: %s\n"%version_n)
				print(version)
				break
			elif arg[i] == "-h" or arg[i] == "--help":
				inst_only = True
				print("Welcome to Analysis plot %s:\n"%version_n)
				print("Copyright (C) 2021  Braga, B. C.\nThis program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; use option '-v' for details.\n")
				print("\nUsage:\n\t-v or --version\t\tprints current version, the python libraries needed and the data format expected.\n")
				print("\t-h or --help\t\tprints this message.\n")
				print("\t-type\tQuantity of files used, for data comparison.\n\t\tone: Normal plot with one data file.\n\t\ttwo: Plots two data files with the same axis information.\n\t\tfour: Plots four data files with the same axis information.\n\t\tallin_one: Plots every dataset given in one XY frame.\n")
				print("\t-eng\tSets default texts language to english. If this flag is not called the texts will be set on portuguese.\n")
				print("\t-fontsize\tInteger value for the fontsize of labels and inplot texts (Default=%d).\n"%fontsize)
				print("\t-index3d\t(Valid only for type 'allin_one') Creates a 3D plot with a index axis separating your datasets.\n")
				print("\t-jpeg\t\t(Recommended when dpi is too big) Saves picture directly to 'Figure.jpeg' instead of generating the interface.\n")
				print("\t-mean\t\t(Valid only for type 'allin_one') Print a horizontal line with the mean value of all the data given.\n")
				print("\t-anatp\t\tAnalysis type:\n\t\trmsd;\n\t\trmsf;\n\t\tradgyr;\n\t\tdist;\n\t\tedcomp (Valid only for types 'one' and 'two') Energy decomposition on the format of 'decode_mmpbsa.py'\n\t\tprot_state (Tested only with type one).\n")
				print("\t-multana\t\t(Has priority over 'anatp', and it's empty by default) Sets one analysis type for each file on the input flag '-i'. The first analysis will be set for the first input file, the second analysis for the second file and so on.\n")
				print("\t-edcut\t\tSets the higher energy limit for highlighting residues on the energy decomposition plot. Default:%s\n"%str(mmpbsa_cut))
				print("\t-edinsetXYpos\t\tSets the inset position relative to the original plot (values between 0 and 1). Ex for 10%% Res axis and 5%% energy axis: -edinsetXYpos 0.1,0.05\n")
				print("\t-edinsetX\t\tSets the residue range for the inset plot. Ex: -edinsetX 169-210\n")
				print("\t-edXtick\t\tSets the residue Ticks size. Ex: -edXtick 10\n")
				print("\t-i\t\tinput data file(s) with a name for the plot (separated by space).\n\t\t\tEg.: -i ph7.00_rmsd.dat pH=7.00\n")
				print("\t-stitle\t\t(Valid only for type 'four') Title for comparison plot.\n")
				print("\t-lblcrd\t(Valid only for type 'four') Label coords.\n")
				print("\t-mlbpos\t(Valid only for type 'allin_one') Label position based on axis percentage. Eg.: -mlbpos 0.5 0.8\n")
				print("\t-fram2time\t(For anatp = rmsd or radgyr) Sets X-axis will change from frame to time (pico seconds). Must inform frame-step conversion.\n\t\t\tEg.: -fram2time 10000\n")
				print("\t-nanosec\t\t(For anatp = rmsd or radgyr) Sets time intervals on X-axis to nanoseconds.\n")
				print("\t-dpi\t\tSets the plot quality (recommended to use only with type=one). Default: 300\n")
				print("\nExamples:\n\t$ python3 Analysis_Plots.py -type one -anatp rmsd -i ph7.00_rmsd.dat Production pH=7.00 -fram2time 10000 -dpi 120\n")
				print("\n\t$ python3 Analysis_Plots.py -type two -anatp rmsf -i ph7.00_rmsf.dat Production pH=7.00 ph8.00_rmsf.dat Production pH=8.00\n")
				print("\n\t$ python3 Analysis_Plots.py -type four -stitle Production -anatp radgyr -lblcrd (16.61,200) -i ph7.00_radgyr.dat pH=7.00 ph8.00_radgyr.dat pH=8.00 ph9.00_radgyr.dat pH=9.00 ph10.00_radgyr.dat pH=10.00 -fram2time 10000 -nanosec\n")
				break

			elif arg[i] == "-debug":
				debug = True
				inst_only = True
				# A flag alone cannot set the loop to jump the next argument 
				#continue
			elif arg[i] == "-type":
				tpe = arg[i+1]
				debug_dic['tpe'] = arg[i+1] 
				continue
			elif arg[i] == "-fontsize":
				fontsize = int(arg[i+1])
				debug_dic['fontsize'] = arg[i+1]
				continue
			elif arg[i] == "-mean":
				mean_flag = True
				debug_dic['mean'] = True
				#continue
			elif arg[i] == "-edcut":
				mmpbsa_cut = int(arg[i+1])
				debug_dic['edcomp cut'] = arg[i+1]
				continue
			elif arg[i] == "-edXtick":
				inset_tick = int(arg[i+1])
				debug_dic['edXtick'] = arg[i+1]
				continue
			elif arg[i] == "-edinsetX":
				#169-210
				range_temp         = arg[i+1].split('-')
				mmpbsa_inset_range = (int(range_temp[0]),int(range_temp[1]))
				debug_dic['edinsetX'] = mmpbsa_inset_range
				continue
			elif arg[i] == "-edinsetXYpos":
				#0.1,0.05
				rangep_temp         = arg[i+1].split(',')
				mmpbsa_inset_XY = (float(range_temp[0]),float(range_temp[1]))
				debug_dic['edinsetXYpos'] = mmpbsa_inset_XY
				continue
			elif arg[i] == "-eng":
				eng = True
				debug_dic['eng'] = True
				#continue
			elif arg[i] == "-index3d":
				tpe = 'allin_one_f3d'
				debug_dic['index3d'] = (True,'tpe:'+tpe)
				#continue
			elif arg[i] == "-jpeg":
				debug_dic["-jpeg"] = (True,file_name)
				if arg[i+1] not in flags:
					file_name = arg[i+1]
					debug_dic["-jpeg"] = (True,file_name)
					continue
				#continue
			elif arg[i] == "-anatp":
				anatp = arg[i+1]
				debug_dic['anatp'] = anatp
				if anatp == 'edcomp':
					mmpbsa = True
					debug_dic['edcomp'] = True		
				continue
			elif arg[i] == "-multana":
				cc = i+1
				if "," in arg[cc] and arg[cc+1] in flags:
					mult_ana_plot = arg[cc].split(",")
					i = cc
				else:
					while cc < len(arg):
						if arg[cc] not in flags:
							mult_ana_plot.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				debug_dic['multana'] = mult_ana_plot	
				continue
			elif arg[i] == "-i":
				default_mod = False
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
							break
						cc_t += 1
					cc = cc_t -1
					if debug:
						print((a_t,b_t))
					if len(b_t) > 0:
						File.append( (a_t, b_t) )
					else:
						i = cc  
						break
					cc += 1
				debug_dic['i'] = File	
				continue
			elif arg[i] == "-stitle":
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
				debug_dic['stitle'] = sup_t
				continue
			elif arg[i] == "-lblcrd":
				#temp_ar = arg[i+1][1:-1].split(',')
				#tempx, tempy = arg
				labx = float(arg[i+1])
				laby = float(arg[i+2])
				debug_dic['lblcrd'] = (labx,laby)
				continue
			elif arg[i] == "-mlbpos":
				multi_label_loc = (float(arg[i+1]),float(arg[i+2]))
				debug_dic['mlbpos'] = multi_label_loc
				continue
			elif arg[i] == "-fram2time":
				fram2time = True
				framstp = float(arg[i+1])
				debug_dic['fram2time'] = (True,framstp)
				continue
			elif arg[i] == "-nanosec":
				nano = True
				debug_dic['nanosec'] = True
				#continue
			elif arg[i] == "-dpi":
				dpi = int(arg[i+1])
				debug_dic['dpi'] = dpi
				continue
			cut +=1

		else: #cut!= i means that the current arg[i] was used in the previous iteration
			cut = i+1

	if debug:
		print("Arguments:",arg)
		print("Input read:",debug_dic)
				
	if default_mod and temp_mod != -1:
		tpe,anatp,File,File2,merge_legend,supertitle,mean_flag,rareplot,multi_label_loc,forced_3Dz,eng,fontsize,forced_mean,fmv,mmpbsa,mmpbsa_inset,mmpbsa_inset_range,mmpbsa_cut,interface = temp_mod
	
	if not inst_only and not version_only:
		if File == []:
			print("No argument given on the flag '-i'\n")
		elif not rareplot:
			#print('files:',File)
			ob4 = Analysis_plot(type=tpe, plot_interface=interface, file_name=file_name, debug=debug_inside, mmpbsa=mmpbsa, mmpbsa_inset=mmpbsa_inset, mmpbsa_inset_XY=mmpbsa_inset_XY, mmpbsa_inset_range=mmpbsa_inset_range, mmpbsa_inset_tick=inset_tick, mmpbsa_cut=mmpbsa_cut, halving_ids=halving_ids, fontsize=fontsize, eng=eng, print_mean=mean_flag, names=File, mult_ana_plot=mult_ana_plot, analysisType=anatp, suptitle=supertitle, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi, label_color=color, merge_legend=merge_legend, multi_merge_label_loc=multi_label_loc, forced_3Dzaxis=forced_3Dz, forced_mean=forced_mean, fmv=fmv)
		else:
			# creating an empty object
			ob4   = Analysis_plot(type='tpe', fontsize=fontsize, eng=eng, names=File2, analysisType=anatp, suptitle=supertitle, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi, label_color=color, merge_legend=merge_legend, multi_merge_label_loc=multi_label_loc, forced_mean=forced_mean, fmv=fmv)
			ob4.X = []
			ob4.Y = []
			
			
			ob4.restriction_break = True
			CpHMDXYlabels         = ob4.XY(files=File2)
			labels_cphmd   = [CpHMDXYlabels[i][1] for i in range(len(CpHMDXYlabels))]
			CpHMDeixox     = CpHMDXYlabels[0][0]
			Xcphmd         = ob4.X
			Ycphmd         = ob4.Y
			
			ob4.X          = []
			ob4.Y          = []
			ob4.frameStep = int(framstp/2)
			if len(File) == 0:
				MDXYlabels = 0
				MDeixox    = CpHMDeixox
				labels_md  = []
			else:
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
				ob4.multi_in_one(labels=dual_labels, label_loc=multi_label_loc, titulo =supertitle, eixoy=ob4.ana_type, eixox=dual_eixox, mean=mean_flag)
			else:
				print("Sheesh!")
