'''Plot: affinity energy 'ae' VS distance 'r' between ser160 hydroxyl and a BHET carbonyl'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)

# Quick sort for the first coodinate of my tupples
def swap(t=[(0,'1'),(-6,'2')], ith = 0, jth = 1):
	if ith != jth:
		temp   = t[ith]
		t[ith] = t[jth]
		t[jth] = temp

def partition(t=[(0,'1'),(-6,'2')], begin= 0, end = 1):
	pivot = t[end]
	rpos  = begin -1
	for i in range(begin,end):
		if t[i][0] < pivot[0]:
			rpos += 1
			swap(t,rpos,i)
	swap(t,rpos+1,end) 
	return 	rpos+1

def SortZerothTupple(t=[(0,'1'),(-6,'2')], begin= 0, end = 1):
	if begin < end:
		pivot = partition(t=t,begin=begin,end=end)

		a = SortZerothTupple(t=t,begin=begin,end=pivot-1)
		b = SortZerothTupple(t=t,begin=pivot+1,end=end)

# Distance betwwen two atoms

def dist(atm1=(1,-2,3), atm2=(0,5,-6)):
	x1, y1, z1 = atm1
	x2, y2, z2 = atm2
	r = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
	return round(r,2)

# List of tupples

def list_AEvsR(ser_hydroxyl = (-20.465,-10.036,7.369), file = 'all_Nice_dock.pdb'):
	'''Creates a filtered pdb and returns a list of affinity vs distance between ser160 hydroxyl and a BHET carbonyl '''
	f = open(file,'r')
	temp = f.readlines()
	f.close()

	models = []
	affinities = []
	rcarbs = [] # 2 per model
	hetatm = []
	f = open("Filtered_%s"%file,'w')

	for i in range(len(temp)):
		if 'model' not in temp[i].lower():
			f.write(temp[i])
			if 'REMARK VINA RESULT:' in temp[i]:
				affinities.append(float(temp[i][20:].split()[0]))
			elif 'HETATM' in temp[i]:
				hetatm.append(temp[i])
			elif 'CONECT' in temp[i]:
				temp_c = temp[i].split()
				if len(temp_c) == 5:
					main_atm = int(temp_c[1])
					main_info = hetatm[main_atm - 1].split()
					alt1 = int(temp_c[2])
					alt1_info = hetatm[alt1 - 1].split()
					alt2 = int(temp_c[3])
					alt2_info = hetatm[alt2 - 1].split()
					alt3 = int(temp_c[4])
					alt3_info = hetatm[alt3 - 1].split()
					c_count = 0
					o_count = 0
					for t in [alt1_info[2],alt2_info[2],alt3_info[2]]:
						if t == 'C':
							c_count += 1
						elif t == 'O':
							o_count += 1
					if main_info[2] == 'C' and c_count == 1 and o_count == 2:
						# main == C carbonyl
						main_xyz = (float(main_info[5]),float(main_info[6]),float(main_info[7]))
						rcarbs.append(dist(atm1=ser_hydroxyl,atm2=main_xyz))
		elif 'model' in temp[i].lower() and 'model' in temp[i+1].lower():
			hetatm = []
			models.append(temp[i]) 
			f.write(temp[i])

	f.close()

	label_out = []
	exc = 0
	h = 0
	for i in range(len(rcarbs)):
		exc +=1
		label_out.append( (models[h],(affinities[h],rcarbs[i])) )
		if exc == 2:
			exc = 0
			h += 1

	return label_out

def plot_AEvsR(dpi=300,file_name='Figure_3.jpeg',pointsize=4,Xg=[],Yg=[],g_color='royalblue',Xw=[],Yw=[],w_color='deepskyblue',Xg_ins=[],Yg_ins=[],mark_color='black',titulo ="",eixoy="Vina score (kcal/mol)",eixox="$r_{hc}\,(\AA)$"):
	#                                      'mediumseagreen'                'deepskyblue'
	fonte=[12,14,16][2]
	fig, ax1 = plt.subplots(dpi=dpi)
	ax1.plot(Xg,Yg,'o',color=g_color,ms=pointsize)
	ax1.plot(Xw,Yw,'o',color=w_color,ms=pointsize)
	ax1.set_xlabel(eixox,fontsize=fonte)
	ax1.set_ylabel(eixoy,fontsize=fonte)
	#ax1.set_xlim(2.5,17)
	ax1.set_ylim(-6.0,-1.8)

	# Create a set of inset Axes: these should fill the bounding box allocated to
	# them.
	ax2 = plt.axes([0,0,1,1])
	# Manually set the position and relative size of the inset axes within ax1
	ip = InsetPosition(ax1, [0.5,0.525,0.5,0.475]) # parant axes; inset X edgde, Y edge, width% to parent, height% to parent
	ax2.set_axes_locator(ip)
	# Mark the region corresponding to the inset axes on ax1 and draw lines
	# linking the two axes.
	# high zorder makes it being drawn above other objects, small zorder are drawn first
	mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec=mark_color,zorder=5)

	ax2.plot(Xg_ins,Yg_ins,'o',color=g_color,ms=pointsize)
	plt.savefig(file_name,bbox_inches='tight')
	#plt.show()

def Main(ref = (-1.062,-12.604,-11.209), arg = ['AvR_DockingPlot.py', 'PDBFILE.pdb'], inset_factor=10):

	# ref := hydroxyl on serine 160 from 6eqe.pdb

	# list of models == [('model i', (float('ae'), float('r')) ), ('model i+1', (float('ae'), float('r')) )]
	data = list_AEvsR(ser_hydroxyl=ref,file=arg[1])
	# Good models
	r_g      = [] # X
	ae_g     = [] # Y

	# Worse models
	r_w      = [] # X
	ae_w     = [] # Y
	bw_limit = -5.5
		
	# Lists for top 5 smallest dist. 'r' and smallest sum of 'r' and 'ae'
	rVmodel    = []
	rPaeVmodel = []

	for i in data:
		rVmodel.append( (i[1][1], i[0], "Ae = %f"%i[1][0]) ) 
		rPaeVmodel.append( (i[1][0]+i[1][1], i[0], "Ae = %f; r_hc = %f"%(i[1][0],i[1][1]), i[1]) )
		if i[1][0] <= bw_limit:
			ae_g.append(i[1][0])
			r_g.append(i[1][1])
		else:
			ae_w.append(i[1][0])
			r_w.append(i[1][1])

	# Top 5 with label
	top        = []
	top_r      = []
		
	print("Finding best 10 fits...\n")
		
	# Quick sorting
	SortZerothTupple(t=rVmodel,begin=0,end=len(rVmodel)-1)
	SortZerothTupple(t=rPaeVmodel,begin=0,end=len(rPaeVmodel)-1)

	Xg_ins = []
	Yg_ins = []
	for i in range(int(len(rPaeVmodel)/inset_factor)):
		if rPaeVmodel[i][3][0] <= bw_limit:
			Yg_ins.append(rPaeVmodel[i][3][0])
			Xg_ins.append(rPaeVmodel[i][3][1])

	if len(rVmodel) < 5:
		top_r = rVmodel
	else:
		top_r = rVmodel[0:5]

	if len(rPaeVmodel) < 5:
		top = rPaeVmodel
	else:
		top = rPaeVmodel[0:5]

	print('\nSmallest 5 "r_hc" found:\n')
	for i in top_r:
		print(i)

	print('\nSmallest 5 sum of absolute values "r_hc plus affinity" found:\n')
	for i in top:
		print(i)
		
	f = open('Best10Fits.txt','w')
	f.write('Smallest 5 "r_hc" found:\n')
	f.write(str(top_r))
	f.write('\n')
	f.write('Smallest 5 sum of absolute values "r_hc plus affinity" found:\n')
	f.write(str(top))
	f.write('\n')
	f.close()

	plot_AEvsR(Xg=r_g,Yg=ae_g,Xw=r_w,Yw=ae_w,Xg_ins=Xg_ins,Yg_ins=Yg_ins)
	
def Help(arg = ['AvR_DockingPlot.py', '-h']):
	print("Please run this code as follows:\n(Option 1): $ pythonCompilerV3+ %s PDBFILE.pdb\n"%arg[0])
	print("(Option 2): $ pythonCompilerV3+ %s PDBFILE.pdb ref Xcoord Ycoord Zcoord"%arg[0])

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv
	
	try:
		if len(arg) == 2:
			Main(arg=arg)
		elif arg[2].lower() == 'ref' and len(arg) == 6:
			Main(ref=(float(arg[3]),float(arg[4]),float(arg[5])),arg=arg)
		else:
			Help()
	except:
		Help()
