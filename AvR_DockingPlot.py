'''Plot: afinity energy 'ae' VS distance 'r' between ser160 hydroxyl and a BHET carbonyl'''

import numpy as np
import matplotlib.pyplot as plt

def dist(atm1=(1,-2,3), atm2=(0,5,-6)):
	x1, y1, z1 = atm1
	x2, y2, z2 = atm2
	r = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
	return round(r,2)

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

def plot_AEvsR(dpi=100,Xg=[],Xw=[],Yg=[],Yw=[],titulo ="Amostras_Docking",eixoy="Afinidade (kcal/mol)",eixox="r_hc (angstrom)"):
    plt.figure(dpi=dpi)
    plt.plot(Xg,Yg,'bo',Xw,Yw,'ro')
    plt.title(titulo)
    plt.ylabel(eixoy)
    plt.xlabel(eixox)
    plt.show()

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv
	
	# list of models == [('model i',('ae','r')), ('model i+1',('ae','r'))]
	if len(arg) == 2:
		data = list_AEvsR(file=arg[1])
	else:
		data = list_AEvsR()

	# Good models
	r_g      = [] # X
	ae_g     = [] # Y
	g_models = []
	# Worse models
	r_w      = [] # X
	ae_w     = [] # Y
	w_models = []

	for i in data:
		if i[1][0] <= -4.5:
			ae_g.append(i[1][0])
			r_g.append(i[1][1])
			g_models.append(i[0])
		else:
			ae_w.append(i[1][0])
			r_w.append(i[1][1])
			w_models.append(i[0])


	plot_AEvsR(Xg=r_g,Yg=ae_g,Xw=r_w,Yw=ae_w)