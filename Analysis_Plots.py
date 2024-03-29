# For xy and xyz plots
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import re
import sys
# For images but not plots
from matplotlib.transforms import IdentityTransform
# For interpolation
from scipy.interpolate import make_interp_spline, BSpline
# Showing the colors
import matplotlib.colors as mcolors
import math
from matplotlib.patches import Rectangle
# Drawing a box on top of the plot
import matplotlib.patches as patches

def plot_colortable(colors, *, ncols=4, sort_colors=True):

	cell_width = 212
	cell_height = 22
	swatch_width = 48
	margin = 12

	# Sort colors by hue, saturation, value and name.
	if sort_colors is True:
		names = sorted( colors, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c))))
	else:
		names = list(colors)

	n = len(names)
	nrows = math.ceil(n / ncols)

	width = cell_width * 4 + 2 * margin
	height = cell_height * nrows + 2 * margin
	dpi = 72

	fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
	fig.subplots_adjust(margin/width, margin/height, (width-margin)/width, (height-margin)/height)
	ax.set_xlim(0, cell_width * 4)
	ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
	ax.yaxis.set_visible(False)
	ax.xaxis.set_visible(False)
	ax.set_axis_off()

	for i, name in enumerate(names):
		row = i % nrows
		col = i // nrows
		y = row * cell_height

		swatch_start_x = cell_width * col
		text_pos_x = cell_width * col + swatch_width + 7

		ax.text(text_pos_x, y, name, fontsize=14,horizontalalignment='left', verticalalignment='center')

		ax.add_patch(Rectangle(xy=(swatch_start_x, y-9), width=swatch_width, height=18, facecolor=colors[name], edgecolor='0.7') )
	plt.savefig("ListOfColors.jpeg",bbox_inches='tight',pad_inches=0.15)
	return "ListOfColors.jpeg"

class Analysis_plot:
	max_num = sys.maxsize
	_Types = ['2ph_2merge','1cmp','four','four_2merge','allin_one']
	def __init__(self, type = 'one', print_mean=False, names=[('file.dat','analysis_title')], analysisType = 'rmsd',mult_ana_plot=[],
	frameToTime=False, frameStep = 5*10**4, timeStep = 0.004, nanosec = False, suptitle='Titulo geral', 
	labelpx = 35.0, labelpy = 0.50, dpi = 100, label_color = 'darkblue', merge_legend = '', multi_merge_label_loc = (0.5,0.5),
	forced_3Dzaxis=['ph','index'][1], forced_mean=False, fmv=4.75, eng=False, fontsize=10, mmpbsa=False,
	mmpbsa_inset=[False,True][0], mmpbsa_inset_range=(169,210), mmpbsa_inset_tick=4, mmpbsa_inset_XY=(0.1,0.125),
	plot_interface=False, mmpbsa_cut=-0.5, halving_ids=[], debug=False, file_name='Figure.jpeg', grid=True,
	vline_color='purple',vline_thickness=0.25,vlines=[132,178,209],bool_shift=True,ID_shift=28,plotid='A', range_freq=(4,5),
	legend_only=False,legend_only_text=[],filter_factor=[20,5.0], ytick_list=[],normal_freq = False,freq_grade=100,interp_degree=3,
	focus_x=False,focus_xij=(30,265),background=False, background_color='black',sharey=False,
	focus_box=True,box_x0=25,box_y0=4,box_xsize=25,box_ysize=2):
		'''Parameters:
		
		mult_ana_plot: Default: Empty ([]). If not empty it will have the same size as 'names' and it will correspond to the analysis on each file in 'names'.

		type: Plot '2ph_2merge' datasets beside each other, 'four'/'four_2merge' datasets in a matrix fashion, multiple datasets in one XY frame ('allin_one') or multiple datasets in different plots but on the same colunm ('1cmp').

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
		self.file_name      = file_name
		self.grid           = grid
		self.max_state      = 0
		self.plotid         = plotid
		self.ytick_list     = ytick_list
		self.filter_yfactor = filter_factor[1]
		self.filter_xfactor = filter_factor[0]
		if bool_shift:
			self.ID_shift = ID_shift
		else:
			self.ID_shift = 0
		self.vlines          = vlines
		self.vline_thickness = vline_thickness
		self.vline_color     = vline_color 
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
		self.plot_interface       = plot_interface
		self.multi_ids            = True
		self.normal_freq          = normal_freq
		self.background           = background
		self.background_color     = background_color
		self.focus_x              = focus_x
		self.focus_xij            = focus_xij
		#interpolation attributes
		self.interp_grade  = 200 # numb. of points generated for the plot
		self.interp_degree = interp_degree # Degree 'k' of the BSpline
		self.range_freq_i  = range_freq[0]
		self.range_freq_f  = range_freq[1]
		#box on the '1cmp' plot
		self.focus_box = focus_box
		self.box_x0    = box_x0
		self.box_y0    = box_y0
		self.box_xsize = box_xsize
		self.box_ysize = box_ysize

		if type in ['1cmp','four']:
			self.multi_ids = False

		if type == 'allin_one' or type == 'allin_one_f3d' or type == '1cmp':
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
		if not legend_only and mult_ana_plot != [] and len(mult_ana_plot) == len(names):
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
						new_list.append(['Decomposição de Energia\n(kcal/mol)','Energy Decomposition\n(kcal/mol)'][self.lang_set])
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
		elif not legend_only:
			if len(mult_ana_plot) >0 and len(names)> 2:
				print("'-multana' doesn't include the analysis for all the files!")
				print("-multana: ",mult_ana_plot)
				print("-i:", names)
				XTlabels = -1
			elif self.mmpbsa:
				XTlabels = self.XY_decomp(files=names)
			else:
				XTlabels = self.XY(files=names,ana_type=self.ana_type)
		else:
			self.legend_frame(text=legend_only_text)

		if self.normal_freq:
			self.norm_freq_label = ['Frequência normalizada (%)','Normalized Frequency'][self.lang_set]
			self.norm_freq(grade=freq_grade)

		if debug:
			print("X:", self.X)
			print("Y:",self.Y)	
		elif not legend_only:
			if XTlabels != -1 and type == '1cmp':
				self.plot_1cmp(X=self.X, Y=self.Y, Xaxis=[XTlabels[i][0] for i in range(len(XTlabels))],
				name= [XTlabels[i][1] for i in range(len(XTlabels))],plts=len(XTlabels),sharey=sharey)
			elif XTlabels != -1 and type == 'four':
				self.plot_four(X=self.X, Y=self.Y, Xaxis=[XTlabels[i][0] for i in range(4)],
				name= [XTlabels[i][1] for i in range(len(XTlabels))])
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
				try:
					mdid = int(re.findall(r'[0-9]+',plot_title)[-1])
				except IndexError:
					if self.multi_ids:
						print("No MD id found on plot title, setting to default.")
					mdid = 1
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
							xname = ['Número do Resíduo','Residue Number'][self.lang_set]
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
							if half_data:
								if x_value%2 ==0:
									continue
								x_value = xcount
								y.append( float(data[1]) )
							else:
								y.append( float(data[1]) )

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
			xname = ['Número do Resíduo','Residue Number'][self.lang_set]
			if not self.mult_ana:
				self.ana_type = ['Decomposição de Energia (kcal/mol)','Energy Decomposition (kcal/mol)'][self.lang_set]
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
					if "Total" not in i and i != '' and i != '\n' and len(data)>=3:
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
				elif abs(false_peak[-1][1]-false_peak[-2][1]) > abs(min(Y))/self.filter_yfactor:
					#i don't want to see a bunch of neighbouring residues 
					notes.append( (X[i],Y[i]) )
				elif abs(false_peak[-1][0]-false_peak[-2][0]) > self.filter_xfactor:
					notes.append( (X[i],Y[i]) )
		return notes	

	def norm_freq(self,grade=100):
		# Apply the normalized frequency for the current analysis (tested only on 'Distance'/Time)  
		# Only applicable(now) if not mult_ana
		X = []
		Y = []
		# self.Y := chosen dataset
		may = max([max(i) for i in self.Y])
		miy = min([min(i) for i in self.Y]) 	
		delta = (may - miy)/grade
		# shared x-axis subdivision
		new_X = [] # new X[i]
		for i in range(int(grade)+1):
			temp = round(miy+i*delta,2)
			new_X.append(temp)
		# counter for each Y dataset on each x region
		freq_count = []
		for i in range(len(self.Y)):
			freq_count.append({})
			for j in new_X:
				freq_count[i][j]=0
		f = open("%s_pars.txt"%self.file_name[:-5],'w')
		for i in range(len(self.Y)):
			# old self.X[i] ::= time (ns)
			# old self.Y[i] ::= distance (Ang)
			data_cnt = 0.0				
			for ii in self.Y[i]:
				#ii: Y-point of the dataset
				for j in freq_count[i]:
					#j: x-region
					if abs(ii-j) <= delta/2:
						data_cnt += 1
						freq_count[i][j]+=1
						break
			freq = [] # new Y[i]
			total_freq = 0
			#range 4-5 ang frequency
			range_freq = 0
			mean_value = 0
			for ii in range(len(new_X)):
				temp = round(freq_count[i][new_X[ii]]/data_cnt,2)
				total_freq += temp
				if new_X[ii] <= self.range_freq_f and new_X[ii] >= self.range_freq_i:
					range_freq += temp
				# para plot em pt-br vou colocar em (%)
				if self.lang_set == 0:
					freq.append( int(temp*100) ) ##
				else:
					freq.append( temp ) ##
				mean_value += new_X[ii]*temp
			# finding max freq x-value
			temp_max = max(freq)
			for c_temp in range(len(freq)):
				if freq[c_temp] == temp_max:
					break
			f.write("Curve %d; Integral: %.2f, freq_%.2f-%.2f: %.2f, Max value: %.2f, Mean: %.2f\n"%(i+1,round(total_freq,2),self.range_freq_i,self.range_freq_f,round(range_freq,2),new_X[c_temp],round(mean_value,2)))
			# Doing interpolation to smooth out the lines
			spl = make_interp_spline(new_X, freq, k=self.interp_degree) #k:= interpolation degree
			X_interp = np.linspace(miy,may,self.interp_grade)
			Y_interp = spl(X_interp)
			for ii in range(len(X_interp)):
				if Y_interp[ii]<0:
					Y_interp[ii] = 0
				'''else:
					# para plot em pt-br vou colocar em (%)
					if self.lang_set == 0:
						Y_interp[ii] = int(Y_interp[ii]*100)'''
			X.append(X_interp)#(new_X)
			Y.append(Y_interp)#(freq)
			
		f.close()
		self.X = X
		self.Y = Y

	def plot_textid(self):
		#finding postion for id
		if self.mult_ana:
			plt.text(0.05, 0.9, '(%s)'%self.plotid, color='steelblue', fontsize=1.5*self.fontsize, transform=plt.gcf().transFigure)
		else:
			max_x, max_y = (-self.max_num, -self.max_num)
			min_x = self.max_num
			if 'list' in str(type(self.Y[0])) or 'numpy' in str(type(self.Y[0])):
				for yy in self.Y:
					max_y = max(max_y,max(yy))
				for xx in self.X:
					max_x = max(max_x,max(xx))
					min_x = min(min_x,min(xx))
			else:
				max_y = max(max_y,max(self.Y))
				max_x = max(max_x,max(self.X))
				min_x = min(min_x,min(self.X))
			id_x = min_x - (max_x - min_x)/5.0
			if 'numpy' in str(type(self.Y[0])):
				id_x += -0.1
			id_y = max_y
			plt.text(x=id_x,y=id_y,s='(%s)'%self.plotid,color='steelblue', fontsize=1.5*self.fontsize)

	def legend_frame(self,text=[r"IQ: $\sigma_i=15$"]):
		fig = plt.figure(dpi=self.dpi)
		cc=0
		'''fig.text(10, (len(text)*100)-cc, '', color='steelblue', fontsize=self.fontsize,
        transform=IdentityTransform())'''
		#500 limit of y
		p = 500 /len(text)
		for i in range(len(text)):
			fig.text(50, 500 -(i+0.5)*p, text[i], color='steelblue', fontsize=self.fontsize,
        	transform=IdentityTransform())
		
		if not self.plot_interface:
			plt.savefig("Legend_box.jpeg",bbox_inches='tight',pad_inches=0.15)
		else:
			plt.show()
	
	def plot_1cmp(self, X = [[], []], Y = [[], []], Xaxis = ["Frames"], name = ["Título"], plts=3, sharey= False):
		'''Plots two X-Y datasets sharing x,y - axis.'''
		
		mark_color  ="lightsteelblue"
		plot_color  ="cornflowerblue"
		inset_color = "darkslategray"
		axs= [0]*plts
		if plts >1:
			fig, axs[0:] = plt.subplots(nrows=plts, ncols=1, sharex=True, sharey= sharey, dpi=self.dpi)
		else:
			fig, ax = plt.subplots(nrows=plts, ncols=1, sharex=True, sharey= sharey, dpi=self.dpi)
			axs = [ax]

		if self.plotid != '':
			self.plot_textid()
		
		if self.background:
			fig.patch.set_facecolor('silver') #'lightsteelblue'

		anatest = ['Estado de Protonação','Protonation State'][self.lang_set]
		anatest2 = ['Decomposição de Energia (kcal/mol)','Energy Decomposition (kcal/mol)'][self.lang_set] 
		if not self.mult_ana:
			if self.ytick_list == []:
				max_y = -self.max_num
				min_y = self.max_num
				for yy in Y:
					max_y = max(max_y,max(yy))
					min_y = min(min_y,min(yy))
				y_tick=np.arange(round(min_y,1),round(max_y,1)+1,round((max_y-min_y)/4,1))
			else:
				y_tick=np.array(self.ytick_list)
				min_y = min(y_tick)
				max_y = max(y_tick)
			plt.setp(axs, yticks=y_tick)
		else:
			jj = -1
			for ts in axs:
				jj+=1 
				if self.ana_list[jj]==anatest:
					y_tick=range(0,self.max_state+1)
				else:
					if self.ytick_list == []:
						max_y = max(Y[jj])
						min_y = min(Y[jj])
						y_tick=np.arange(round(min_y,1),round(max_y,1)+1,round((max_y-min_y)/3,1))
					else:
						y_tick=np.array(self.ytick_list)
						min_y = min(y_tick)
						max_y = max(y_tick) 
				plt.setp(ts, yticks=y_tick)
		plt.subplots_adjust(hspace=0.4)

		diff_y5t = len(y_tick) -5
		data_i = -1
		inset  = []
		bool_erf = self.mmpbsa or 'RMSF' in self.ana_type
		for ts in axs:
			# to draw something (a box in this case) on the plot
			if self.focus_box:
				box = patches.Rectangle((self.box_x0, self.box_y0), self.box_xsize, self.box_ysize, linewidth=1, edgecolor='r', facecolor='none',zorder=len(axs)+1)
				ts.add_patch(box)
			if self.background:
				ts.set_facecolor(self.background_color)
			if self.focus_x:
				ts.set_xlim(self.focus_xij[0], self.focus_xij[1]) #(30,265)

			if diff_y5t>0:
				ts.tick_params(axis='y', labelsize=(1.0-0.1*diff_y5t)*self.fontsize)
			data_i +=1
			bool_me  = self.mult_ana and self.ana_list[data_i]=='edcomp'
			bool_mrf = self.mult_ana and 'rmsf' in self.ana_list[data_i]
			if bool_me or bool_mrf or bool_erf:
				X_0 = []
				for xx in X[data_i]:
					X_0.append(xx+self.ID_shift)
			else:
				X_0 = X[data_i]

			if self.mult_ana and self.ana_list[data_i]==anatest or anatest in self.ana_type:
				ts.plot(X_0,Y[data_i],'o',ms=1)
			elif anatest2 in self.ana_type:
				print('Residue correction for plot %s:'%(data_i+1), self.ID_shift)
				ts.plot(X_0,Y[data_i],color=plot_color)
				ts.set_ylim(round(min_y,1)-0.5,round(max_y,1)+1)
			else:
				ts.plot(X_0,Y[data_i],color=plot_color)
			ts.set_title(name[data_i], fontsize=self.fontsize)
			ts.grid(self.grid,color='gray',linewidth=0.5)
			
			#esboco zuado
			##ylb='Binding Energy' ##########
			##plt.setp(ts, yticks=[0],xticks=[])#######
			##ts.text(160, -2, 'S', color='orange', fontsize=self.fontsize)
			#ts.text(161, -2.5, 'M', color='blue', fontsize=self.fontsize)
			#if ylb=='Binding Energy' and ts == axs[1]:
			#	ts.text(237, -1, 'K', color='green', fontsize=self.fontsize)

			if self.mult_ana and self.ana_list[data_i]=='edcomp' or self.mmpbsa:
				# inset data res Eg:170-210
				if self.mmpbsa_inset:
					inset1_x   = []
					inset1_y   = []
					for i in self.mmpbsa_inset_range:
						inset1_x.append( X_0[i] )
						inset1_y.append( Y[data_i][i] )
					temp = plt.axes([0,0,1,1],label='upinset')
					inset.append(temp)	
					inset[data_i].tick_params(axis='both', which='major', labelsize=0.5*self.fontsize)
					inset[data_i].set_xticks(np.arange(self.mmpbsa_inset_range[0], self.mmpbsa_inset_range[-1]+1, self.mmpbsa_inset_restick), minor=False)
					ip = InsetPosition(ts, [self.mmpbsa_inset_XY[0],self.mmpbsa_inset_XY[1],0.4,0.45])
					inset[data_i].set_axes_locator(ip)
					mark_inset(ts, inset[data_i], loc1=2, loc2=4, fc="none", ec=mark_color,zorder=-1)
					inset[data_i].plot(inset1_x,inset1_y,color=plot_color)
					#inset1.axvline(x=178,color='green')
					#inset1.axvline(x=209,color='red')
			
				note = self.peaks(X=X_0,Y=Y[data_i])
				for i in note:
					def_weight='normal'
					if self.background and self.background_color=='black':
						text_color  ='white'
						self.vline_color = 'gold'
					else:
						text_color  ='black'
					cor   = 0
					if i[1] >0:
						if self.background and self.background_color=='black':
							text_color = 'orange'#'yellow'
							def_weight = 'bold'
						else:
							text_color = 'darkred'
						cor   = [0,0.2,0.4][2]
					if i[0] == 160:
						cor = 0.5
					elif i[0] == 206:
						cor = -0.5
					'''# do esboço zuado!!
					if ylb!='Binding Energy':
						nota = self.mmpbsa_res[i[0]-(1+self.ID_shift)]+str(i[0])
					elif i[0] in res_petase:
						text_color = res_petase[i[0]][0]
						nota = res_petase[i[0]][1]
					elif i[0]==206:
						if ts == ax1:
							text_color = 'orange'
							nota = 'D'
						else:
							text_color = 'green'
							nota = 'E'
					elif i[0]==237:
						if ts == ax1:
							text_color = 'red'
							nota = 'H (Unfavorable)' '''
					if self.mmpbsa_inset:
						if i[0] in self.mmpbsa_inset_range:
							inset[data_i].text(i[0],i[1],self.mmpbsa_res[i[0]-(1+self.ID_shift)]+str(i[0]), fontweight=def_weight, color=inset_color, fontsize=self.fontsize*0.5)
						else:
							ts.text(i[0], i[1]-cor, self.mmpbsa_res[0][i[0]-(1+self.ID_shift)]+str(i[0]), fontweight=def_weight, color=text_color, fontsize=self.fontsize*0.75)
					else:		
						ts.text(i[0], i[1]-cor, self.mmpbsa_res[data_i][i[0]-(1+self.ID_shift)]+str(i[0]), fontweight=def_weight, color=text_color, fontsize=self.fontsize*0.75)

		if self.mult_ana:
			counter = -1
			for ts in axs:
				counter += 1
				ts.set_xlabel(Xaxis[counter],fontsize=self.fontsize)
				ts.set_ylabel(self.ana_list[counter],fontsize=self.fontsize*0.8)	
		elif self.mmpbsa:
			counter=-1			
			for ts in axs:
				counter+=1
				ts.set_xlabel(Xaxis[counter],fontsize=self.fontsize)
				for xid in self.vlines:
					ts.axvline(x=xid+self.ID_shift,color=self.vline_color,lw=self.vline_thickness,zorder=-1)
			if plts==1:
				axs[0].set_ylabel(self.ana_type,fontsize=self.fontsize)
			elif plts==2:
				axs[1].set_ylabel(self.ana_type,fontsize=self.fontsize*1.25)
				axs[1].yaxis.set_label_coords(-0.07, 1.25)
			elif plts==3:
				axs[1].set_ylabel(self.ana_type,fontsize=self.fontsize*1.25)
			elif plts>3:
				axs[2].set_ylabel(self.ana_type,fontsize=self.fontsize*1.5)
				axs[2].yaxis.set_label_coords(-0.07, 1.25) #1 right above the third plot
		else:
			counter = -1
			for ts in axs:
				counter += 1
				ts.set_xlabel(Xaxis[counter], fontsize=self.fontsize)
				if plts < 3:
					ts.set_ylabel(self.ana_type, fontsize=self.fontsize)
			if plts==3:
				axs[1].set_ylabel(self.ana_type,fontsize=self.fontsize*1.25)
			elif plts>3:
				axs[2].set_ylabel(self.ana_type,fontsize=self.fontsize*1.5)
				axs[2].yaxis.set_label_coords(-0.1, 1)
	
		#pyplot.subplots can hide redundant axes
		for ts in axs:
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
		if self.plotid != '':
			self.plot_textid()
		big_pp = 0
		axs = [ax1, ax2]
		for ap in axs:
			##
			ap.set_ylim(int(min(Y[big_pp]))-2,int(max(Y[big_pp])+2))
			plt.setp(ap, yticks=range( int(min(Y[big_pp]))-1 , int(max(Y[big_pp]))+2))
			##
			ap.plot(X[big_pp],Y[big_pp], label= name[big_pp]) #self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= name[big_pp+1])
			ap.set_xlabel(Xaxis,fontsize=self.fontsize)
			ap.set_ylabel(self.ana_type,fontsize=self.fontsize)
			ap.legend(fontsize=self.fontsize)
			ap.set_title("pH %d"%(big_pp+7), fontsize=self.fontsize)
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid(self.grid)
			for xid in self.vlines:
				ap.axvline(x=xid+self.ID_shift,color=self.vline_color,lw=self.vline_thickness,zorder=-1)
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
		if self.plotid != '':
			self.plot_textid()
		ax1 = plt.subplot(221)
		plt.plot(X[0],Y[0])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[0], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[0], fontsize=self.fontsize)
		plt.grid(self.grid)
		ax2 = plt.subplot(222, sharey=ax1)		
		plt.plot(X[1],Y[1])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[1], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[1], fontsize=self.fontsize)
		plt.grid(self.grid)
		ax3 = plt.subplot(223, sharex=ax1)
		plt.plot(X[2],Y[2])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[2], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[2], fontsize=self.fontsize)
		plt.grid(self.grid)
		ax4 = plt.subplot(224, sharex=ax2, sharey=ax3)
		plt.plot(X[3],Y[3])
		plt.text(x=self.labelpx, y=self.labelpy, s=name[3], color=self.label_color, fontsize=self.fontsize)
		plt.ylabel(self.ana_type, fontsize=self.fontsize)
		plt.xlabel(Xaxis[3], fontsize=self.fontsize)
		plt.grid(self.grid)
		if self.supertitle:
			plt.suptitle(self.suptitle,fontsize=self.fontsize)
		
		for ts in [ax1,ax2,ax3,ax4]:
			ts.label_outer()
			for xid in self.vlines:
				ts.axvline(x=xid+self.ID_shift,color=self.vline_color,lw=self.vline_thickness,zorder=-1)

		if not self.plot_interface:
			plt.savefig(self.file_name,bbox_inches='tight')
		else:
			plt.show()

	def plot_four_2merge(self, X = [[], [], [], []], Y = [[], [], [], []], Xaxis = "Frames", name = ["Título"]):
		'''Plots two X-Y datasets sharing x,y - axis.'''

		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, dpi=self.dpi)
		if self.plotid != '':
			self.plot_textid()
		big_pp = 0
		axs = [ax1, ax2, ax3, ax4]
		for ap in axs:
			ap.plot(X[big_pp],Y[big_pp], label= self.merge_legend[0])
			ap.plot(X[big_pp+1],Y[big_pp+1], label= self.merge_legend[1])
			ap.legend(fontsize=self.fontsize)
			ap.set_title(name[big_pp], fontsize=self.fontsize)
			#ap.text(x=self.labelpx, y=self.labelpy, s=name[big_pp], color=self.label_color)
			ap.grid(self.grid)
			for xid in self.vlines:
				ap.axvline(x=xid+self.ID_shift,color=self.vline_color,lw=self.vline_thickness,zorder=-1)
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

		if self.plotid != '':
			self.plot_textid()
		
		if not self.background:
			plt.figure(dpi=self.dpi)
		else:
			fig = plt.figure(dpi=self.dpi)
			fig.patch.set_facecolor('silver') #'lightsteelblue'
			ts = fig.add_subplot(1, 1, 1)
			ts.set_facecolor(self.background_color)

		
		plt.grid(self.grid)
		mean_value = []
		for i in range(len(self.X)):
			mean_value.append( sum(self.Y[i])/len(self.Y[i]) )
			if self.mmpbsa or 'RMSF' in self.ana_type:
				X_0 = []
				for xx in self.X[i]:
					X_0.append(xx+self.ID_shift)
			else:
				X_0 = self.X[i]
			#plt.plot(X_0, self.Y[i],'-o',ms=1) # ponto	
			plt.plot(X_0, self.Y[i]) # linha

		for xid in self.vlines:
			plt.axvline(x=xid+self.ID_shift,color=self.vline_color,lw=self.vline_thickness,zorder=-1)
			
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
		
		if not self.normal_freq:
			plt.ylabel(eixoy, fontsize=self.fontsize)
			plt.xlabel(eixox, fontsize=self.fontsize)
		else:
			plt.ylabel(self.norm_freq_label, fontsize=self.fontsize)
			plt.xlabel(self.ana_type, fontsize=self.fontsize)
			
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
		if self.plotid != '':
			self.plot_textid()
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
	arg = sys.argv
	version_n = 1.7
	version = ''' Analysis plot - Pyplot extension for RMSD, RMSF and Radgyr data analysis.

Python3 modules needed: 
	numpy; 
	matplotlib;
	mpl_toolkits; # i believe this comes together with matplotlib
	re;           # i believe this is a default python library but i could be wrong
	scipy.


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
	dic_Type = {13:'allin_one', 133:'allin_one_f3d', -13:'1cmp', 22: '2ph_2merge', 4:'four', 42: 'four_2merge'}
	
	#default keys
	sharey       = False # for type '1cmp' only
	interface    = True
	mmpbsa             = False
	mmpbsa_inset       = False
	mmpbsa_inset_range = (169,210)
	mmpbsa_cut         = -0.5
	version_only = False
	inst_only    = False
	anatp        = ['rmsd','rmsf','radgyr','dist','prot_state'][3]
	tpe          = dic_Type[13]
	fontsize     = 10
	eng          = False
	forced_3Dz   = 'index'
	color        = 'black' #'red' #'darkblue' #label only
	mut          = ['HIS 237 - GLU','ASP 206 - GLU'][1]
	supertitle   = '' #['Produção - Mutação %s'%mut,'Petase nativa - Produção'][1]
	bckbne_comp  = False 
	merge_legend = ''
	labx, laby   = (100.5, 0.51) #rmsd
	fram2time    = False
	framstp      = 62500 # Ultra #gpu-high: 50 000
	nano         = True  # False
	dpi          = [100,200,300][2]
	grid         = True
	mean_flag    = False
	normal_freq  = False
	freq_grade   = 20
	# drawing box
	focus_box = False
	box_x0    = 25
	box_y0    = 4
	box_xsize = 25
	box_ysize = 2

	# for interpolation
	interp_degree = 3
	range_freq    = (4,5)

	# debug
	debug        = False # for the arguments only
	debug_dic    = {}

	# special cases #IGNORE THIS
	File         = []
	File2        = []
	merge_legend = ['NativaS1',['NativaS2','D206E'][1]]
	forced_mean  = False
	fmv          = 4.75
	ytick_list   = []
	rareplot     = False
	multi_label_loc = (.685,0.5)
	#mmpbsa_cut must be negative!!
	mult_ana_plot   = [] # multiple analysis
	vline_color     = 'purple'
	vline_thickness = 0.25
	vlines          = [] #[132,178,209]
	filter_factor   = [20,5.0]
	bool_shift      = False #ver flag -xshift
	ID_shift        = 28
	plotid          = ''
	file_name       = 'Figure.jpeg'
	focus_x          = False
	focus_xij        = (30,265)
	background       = False
	background_color = 'black'
	# for legend box only and no plots
	## on flag '-i' i've set this to False so it does
	#  not matter even if it's True here
	l_o   = True 
	# text to be printed
	l_o_t = ['(A) Nativa','(B) H237K','(C) D206E/H237K']
	#['(A) Wildtype','(B) H237K','(C) D206E/H237K']
	#['(A) All','(B) Wildtype','(C) D206E','(D) D206E/H237K','(E) H237K']

	default_mod     = False
	inset_tick      = 10
	mmpbsa_inset_XY = (0.10,0.10)
	halving_ids = [[],[1,2,3]][0]
	met= 0 # MD ## met=1 := CpHMD
	# If the flag '-i' is read then 'default_mod' becames False  
	#tpe=-13 eh 1 coluna 3 linhas
	temp_mod = Default_modifier(on=default_mod,tpe=dic_Type[13],
	anatp=['rmsd','rmsf','radgyr','dist','prot_state'][3],File=File,File2=File2, 
	enzyme=['Nat','D206E','D206EH237K','H237K','H237Knored','Nat_nored','WT-D-DH'][0], mmpbsa=[True,False][1],
	met=0,mmpbsa_cut=[-0.5,-1.5][0], mmpbsa_inset=[True,False][1], mmpbsa_inset_range=[(169,210),(150,160)][0],
	merge_legend=merge_legend,supertitle=['%s-BHET pH7'%['IsPETase','D206EH237K'][1],''][1],
	mean_flag=[True,False][1],rareplot=[True,False][1],interface=[True,False][1],
	multi_label_loc=(.685,0.5),forced_3Dz=['ph','index'][0],
	forced_mean=[True,False][1], fmv=[4.65,4.75,4.97][0], ph_reverse=[True,False][1],
	eng=[True,False][0],fontsize=[8,12,14][1],mdlistrange=range(1,6),
	igph7md=[],igph9md=[1,2,3,4,5,6],igph7cphmd=[],igph9cphmd=[],
	res=[('GL4178','GLU 206'),('LYS209','LYS 237'),('AS4178','ASP 206'),('HIP209','HIS 237'),('TYR59','TYR 87')][3])

	#print("\nmodifiers:",temp_mod)
	
	if anatp == 'rmsf':
		fram2time  = False
		nano       = False
		labx, laby = (150, 2.25) 
	elif anatp == 'radgyr':
		labx, laby = (130, 16.8)

	flags = ["&","-v","--version","-redbox","-backgrnd","-yshare","-xlimited","-normfreq","-rangfreq","-splinedegree","-ytick","-id","-ffa","-ffax","-hidden","-xshift","-vcolor","-vlines","-vltcs","-fontsize","-jpeg","-nogrid","-edinsetXYpos","-edinsetX","-edXtick","-edcut","-eng","-h","--help","-type","-index3d","-mean", "-i", "-debug","-anatp","-multana","-stitle","-fram2time", "-framstp","-nanosec","-dpi","-lblcrd","-mlbpos"]

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
				print("\t-type\tQuantity of files used, for data comparison.\n\t\t1cmp: 1 column and multiple plots (Eg.: '-type 1cmp -i 1.dat name1' for one and '-type 1cmp -i 1.dat name1 2.dat name2 3.dat name3' for three plots in one column).\n\t\tfour: Plots four data files with the same axis information.\n\t\tallin_one: Plots every dataset given in one XY frame.")
				print("\t-eng\tSets default texts language to english. If this flag is not called the texts will be set on portuguese.\n")
				print("\t-vlines\tPlot vertical lines. Eg.: -vlines 178,209\n")
				print("\t-vcolor\tSets the color for the vertical lines (Default: purple).\n")
				print("\t-vltcs\tSets the thickness of the vertical lines (Default: 0.25).\n")
				print("\t-ffa\t(Valid for 'edcomp' only) Filter factor for annotations (Default: 5). A small value means a bigger y-distance between annotations, so it will hide more annotations. Useful to change it when there are many enrgy wells close to one another.\n")
				print("\t-redbox\t(Valid only for the '1cmp' type) Draw a red box on all plots sharing its parameters: initial x (x0), initial y (y0), length of x (xl) and length of y (yl) [MUST be given in this order]\n\t\t\tEg.: -redbox 25 4 23 6.\n")
				print("\t-yshare\t(Valid only for the '1cmp' type) Forces all the plots to share the zoom on the Y-axis.\n")
				if "-hidden" in arg:
					color_file = plot_colortable(mcolors.CSS4_COLORS)
					print("\t-backgrnd\t(Valid only for the '1cmp'  and 2D 'allin_one' type) Changes the color for the background. Eg: \"-backgrnd black\". The list of colors are presented on file %s\n"%color_file)
					print("\t-xlimited\t(Valid only for the '1cmp' type) Focus the X-axis on a specific range. Eg: \"-xlimited 30 265\" could be used for rmsf or edcomp to reduced the residues shown.\n")
					print("\t-normfreq\t(Valid only for the usual 2D style 'allin_one' with 'rmsd','radgyr' or 'dist') Sets analysis as a normalized frequency for the Y-axis. If a number of subdivions is not given then this program will set it to %d by default.\tPs: We will ignore the negative parts of our interpolation function.\n"%freq_grade)
					print("\t-rangfreq\t(Valid only for the usual 2D style 'allin_one' with 'rmsd','radgyr' or 'dist') Sets a range to calculate the frequency. Eg: \"-rangfreq 4 5\" shows the frequency the data stays between 4 and 5 of the X-axis of the normalized frequency plot.\n")
					print("\t-splinedegree\t(Valid only for 'normfreq') Sets the interporlation degree (Default: 3). Must be a natural number.\n")
					print("\t-ffax\t(Valid for 'edcomp' only) X-axis filter factor for annotations (Default: 20). Differently from 'ffa' it means the exact number of residues it will skip if 'ffa' fails.\n")
					print("\t-ytick\t(Valid for '1cmp' only) Y-axis tick list.\n")
				print("\t-xshift\t(Valid for 'edcomp' and 'rmsf' only) Adds to the X points, to shift the data in the x-axis (Default: 0).\n")
				print("\t-fontsize\tInteger value for the fontsize of labels and inplot texts (Default=%d).\n"%fontsize)
				print("\t-index3d\t(Valid only for type 'allin_one') Creates a 3D plot with a index axis separating your datasets.\n")
				print("\t-jpeg\t\t(Recommended when dpi is too big) Saves picture directly to 'Figure.jpeg' instead of generating the interface.\n")
				print("\t-nogrid\t\tTurn off 'grid' option.\n")
				print("\t-mean\t\t(Valid only for type 'allin_one') Print a horizontal line with the mean value of all the data given.\n")
				print("\t-anatp\t\tAnalysis type:\n\t\trmsd;\n\t\trmsf;\n\t\tradgyr;\n\t\tdist;\n\t\tedcomp (Valid only for types 'one' and '1cmp') Energy decomposition on the format of 'decode_mmpbsa.py'\n\t\tprot_state (Tested only with type one).\n")
				print("\t-multana\t\t(Has priority over 'anatp', and it's empty by default) Sets one analysis type for each file on the input flag '-i'. The first analysis will be set for the first input file, the second analysis for the second file and so on.\n")
				print("\t-edcut\t\tSets the higher energy limit for highlighting residues on the energy decomposition plot. Default:%s\n"%str(mmpbsa_cut))
				print("\t-edinsetXYpos\t\tSets the inset position relative to the original plot (values between 0 and 1). Eg. for 10%% Res axis and 5%% energy axis: -edinsetXYpos 0.1,0.05\n")
				print("\t-edinsetX\t\tSets the residue range for the inset plot. Eg.: -edinsetX 169-210\n")
				print("\t-edXtick\t\tSets the residue Ticks size. Ex: -edXtick 10\n")
				print("\t-i\t\tInput data file(s) with a name for the plot (separated by space).\n\t\t\tEg.: -i ph7.00_rmsd.dat pH=7.00\n")
				print("\t-id\t\tAn identification number or letter for the plot\n")
				print("\t-stitle\t\tTitle for comparison plot.\n")
				print("\t-lblcrd\t(Valid only for type 'four') Label coords.\n")
				print("\t-mlbpos\t(Valid only for type 'allin_one') Label position based on axis percentage. Eg.: -mlbpos 0.5 0.8\n")
				print("\t-fram2time\t(For anatp = rmsd or radgyr) Sets X-axis will change from frame to time (pico seconds). Must inform frame-step conversion.\n\t\t\tEg.: -fram2time 10000\n")
				print("\t-nanosec\t\t(For anatp = rmsd or radgyr) Sets time intervals on X-axis to nanoseconds.\n")
				print("\t-dpi\t\tSets the plot quality (recommended to use only with type=one). Default: 300\n")
				print("\nExamples:\n\t$ python3 Analysis_Plots.py -type 1cmp -anatp rmsd -i ph7.00_rmsd.dat Production pH=7.00 -fram2time 10000 -dpi 120\n")
				print("\n\t$ python3 Analysis_Plots.py -type 1cmp -anatp rmsf -i ph7.00_rmsf.dat Production pH=7.00 ph8.00_rmsf.dat Production pH=8.00\n")
				print("\n\t$ python3 Analysis_Plots.py -type four -stitle Production -anatp radgyr -lblcrd (16.61,200) -i ph7.00_radgyr.dat pH=7.00 ph8.00_radgyr.dat pH=8.00 ph9.00_radgyr.dat pH=9.00 ph10.00_radgyr.dat pH=10.00 -fram2time 10000 -nanosec\n")
				break
			# THERE ARE MANY DUPLICATED FLAG-VERIFICATIONS BECAUSE
			# I DIDN'T NEED THIS TO BE OPTMIZED, AND ONLY CARED IF IT WAS WORKING
			# WE CAN MAKE A FEW METHODS TO SIMPLIFY THE IMPLEMENTATION BELOW!
			elif arg[i] == "-debug":
				debug = True
				#inst_only = True
				# A flag alone cannot set the loop to jump the next argument 
				#continue
			elif arg[i] == "-type":
				tpe = arg[i+1]
				debug_dic['tpe'] = tpe 
				continue
			elif arg[i] == "-backgrnd":
				background = True
				background_color = arg[i+1]
				debug_dic['backgrnd'] = (background,background_color) 
				continue
			elif arg[i] == "-xlimited":
				focus_x = True
				#4.1,5
				cc = i+1
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
					temp = arg[cc].split(",")
					i = cc
				else:
					# 4.1 5
					temp = []
					while cc < len(arg):
						if arg[cc] not in flags:
							temp.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				
				if len(temp) !=2:
					inst_only = True
					print("Expected 2 values on this argument! Given:", temp)
					break
				focus_xij = (float(temp[0]),float(temp[1]))
				debug_dic['xlimited'] = (True, focus_xij)		
				continue
			elif arg[i] == "-splinedegree":
				interp_degree = int(arg[i+1])
				debug_dic['splinedegree'] = interp_degree 
				continue
			elif arg[i] == "-normfreq":
				normal_freq = True
				if i+1< len(arg) and arg[i+1] not in flags:
					freq_grade  = float(arg[i+1])
					debug_dic['normfreq'] = (True,freq_grade)
					continue
				else:
					debug_dic['normfreq'] = (True,freq_grade)
			elif arg[i] == "-rangfreq":
				#4.1,5
				cc = i+1
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
					temp = arg[cc].split(",")
					i = cc
				else:
					# 4.1 5
					temp = []
					while cc < len(arg):
						if arg[cc] not in flags:
							temp.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				
				if len(temp) !=2:
					inst_only = True
					print("Expected 2 values on this argument! Given:", temp)
					break
				range_freq = (float(temp[0]),float(temp[1]))
				debug_dic['rangfreq'] = range_freq 
				continue
			elif arg[i] == "-ffa":
				filter_factor[1] = float(arg[i+1])
				debug_dic['ffa'] = filter_factor[1] 
				continue
			elif arg[i] == "-ffax":
				filter_factor[0] = float(arg[i+1])
				debug_dic['ffax'] = filter_factor[0] 
				continue
			elif arg[i] == "-id":
				plotid = arg[i+1]
				debug_dic['id'] = plotid 
				continue
			elif arg[i] == "-xshift":
				ID_shift   = int(arg[i+1])
				bool_shift = True
				debug_dic['xshift'] = ID_shift 
				continue
			elif arg[i] == "-vcolor":
				vline_color = arg[i+1]
				debug_dic['vcolor'] = vline_color 
				continue
			elif arg[i] == "-ytick":
				ntype =[int,float]
				ytty  = 0 #Default is integer
				#2,0,-2,-4,-6
				cc = i+1
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
					if '.' in arg[cc]:
						ytty = 1
					temp = arg[cc].split(",")
					i = cc
				else:
					# 2 0 -2 -4 -6
					temp = []
					while cc < len(arg):
						if arg[cc] not in flags:
							if '.' in arg[cc]:
								ytty = 1
							temp.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				ytick_list = []
				for jj in temp:
					ytick_list.append(ntype[ytty](jj))
				ytick_list.sort()
				debug_dic['ytick'] = ytick_list 
				continue
			elif arg[i] == "-vlines":
				#132,178,209
				cc = i+1
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
					temp = arg[cc].split(",")
					i = cc
				else:
					# 132 178 209
					temp = []
					while cc < len(arg):
						if arg[cc] not in flags:
							temp.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				vlines = []
				for jj in temp:
					vlines.append(int(jj))
				debug_dic['vlines'] = vlines 
				continue
			elif arg[i] == "-vltcs":
				vline_thickness  = float(arg[i+1])
				debug_dic['vltcs'] = vline_thickness
				continue
			elif arg[i] == "-fontsize":
				fontsize = int(arg[i+1])
				debug_dic['fontsize'] = arg[i+1]
				continue
			elif arg[i] == "-mean":
				mean_flag = True
				debug_dic['mean'] = True
				#continue
			elif arg[i] == "-yshare":
				sharey = True
				debug_dic['yshare'] = True
				#continue
			elif arg[i] == "-edcut":
				mmpbsa_cut = float(arg[i+1])
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
			elif arg[i] == "-nogrid":
				grid = False
				debug_dic['grid'] = grid
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
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
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
			elif arg[i] == "-redbox":
				focus_box = True
				cc = i+1
				bool_check = "," in arg[cc]
				if cc+1 < len(arg):
					bool_check = "," in arg[cc] and arg[cc+1] in flags
				
				if bool_check:
					boxpar = arg[cc].split(",")
					i = cc
				else:
					boxpar = []
					while cc < len(arg):
						if arg[cc] not in flags:
							boxpar.append( arg[cc] )
						else:
							i = cc - 1
							break
						cc += 1
				if len(boxpar) != 4:
					inst_only = True
					print("With flag '-redbox', you must give FOUR numbers representing: 'x0'; 'y0'; 'xsize' and 'ysize' in this order!")
					break
				box_x0    = float(boxpar[0])
				box_y0    = float(boxpar[1])
				box_xsize = float(boxpar[2])
				box_ysize = float(boxpar[3])
				debug_dic['redbox'] = boxpar	
				continue
			elif arg[i] == "-i":
				default_mod = False
				l_o         = False
				File        = []
				#File.append( (arg[i+1], arg[i+2]) ) # arg[i+2] != 'Production pH=7.00' !!
				cc = i+1
				while cc < len(arg):
					#a_t is the file name
					if arg[cc] not in flags:
						a_t = arg[cc]
					else:
						i = cc - 1
						break
					
					#b_t is the label given to the file
					b_t = ''
					cc_t = cc +1
					while cc_t < len(arg):
						if arg[cc_t] not in flags:
							if '.dat' not in arg[cc_t]:
								if len(b_t) > 0:
									b_t += ' '
								b_t += arg[cc_t]
							else:
								if cc_t == cc+1:
									#2 files and no label
									print("\tYou must give a label to your first plot!")
								elif '.dat' in arg[cc_t] and '.dat' in arg[cc_t-1]:
									print("\tYou forgot to label your file!") 
								break
						else:
							break
						cc_t += 1
					cc = cc_t -1
					if debug:
						print((a_t, b_t))
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
			print("No 'complete' argument given on the flag '-i'\n")
		elif not rareplot:
			#print('files:',File)
			ob4 = Analysis_plot(focus_box=focus_box,box_x0=box_x0,box_y0=box_y0,box_xsize=box_xsize,box_ysize=box_ysize,sharey=sharey,focus_x=focus_x,focus_xij=focus_xij,background=background,background_color=background_color,range_freq=range_freq,interp_degree=interp_degree,normal_freq=normal_freq,freq_grade=freq_grade,legend_only=l_o,legend_only_text=l_o_t,ytick_list=ytick_list,filter_factor=filter_factor,type=tpe,plotid=plotid,vline_color=vline_color,vline_thickness=vline_thickness,vlines=vlines,bool_shift=bool_shift,ID_shift=ID_shift,plot_interface=interface, grid=grid, file_name=file_name, debug=debug, mmpbsa=mmpbsa, mmpbsa_inset=mmpbsa_inset, mmpbsa_inset_XY=mmpbsa_inset_XY, mmpbsa_inset_range=mmpbsa_inset_range, mmpbsa_inset_tick=inset_tick, mmpbsa_cut=mmpbsa_cut, halving_ids=halving_ids, fontsize=fontsize, eng=eng, print_mean=mean_flag, names=File, mult_ana_plot=mult_ana_plot, analysisType=anatp, suptitle=supertitle, frameToTime=fram2time, frameStep=framstp,  nanosec=nano, labelpx=labx, labelpy=laby, dpi=dpi, label_color=color, merge_legend=merge_legend, multi_merge_label_loc=multi_label_loc, forced_3Dzaxis=forced_3Dz, forced_mean=forced_mean, fmv=fmv)
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
