import random

random.seed()

#'python3 verify_log.py Log_R%d.out'%i
newpy = open('verify_log.py','w')
newpy.write('if __name__ == "__main__":\n')
newpy.write('\timport sys\n')
newpy.write('\tfrom os import system as cmd\n')
newpy.write('\targ = sys.argv\n')
newpy.write('\tif len(arg) == 2 and "log" in arg[1].lower():\n')
newpy.write('\t\ta1 = arg[1].find("_R")+2\n')
newpy.write('\t\ta2 = arg[1].find(".out")\n')
newpy.write('\t\tinp_file = "config4vina_%s.inp"%(arg[1][a1:a2])\n')###
newpy.write('\t\tg = open(inp_file,"r")\n')
newpy.write('\t\tinp_lines = g.readlines()\n')
newpy.write('\t\tg.close()\n')
newpy.write('\t\tfor temp in inp_lines:\n')
newpy.write('\t\t\tif "out" in temp:\n')
newpy.write('\t\t\t\tdeez = temp.split()\n')
newpy.write('\t\t\t\tout_file = deez[len(deez)-1]\n')###
newpy.write('\t\t\t\tbreak\n')
newpy.write('\t\tf = open(arg[1],"r")\n')
newpy.write('\t\tflist = f.readlines()\n')
newpy.write('\t\tf.close()\n')
newpy.write('\t\tfor i in flist:\n')
newpy.write('\t\t\tdata = i.split()\n')
newpy.write('\t\t\tif len(data) == 4 and "1" in data[0]:\n')
newpy.write('\t\t\t\tbreak\n')
newpy.write('\t\tif float(data[1]) >=0:\n') # affinity >= 0 kcal/mol :: bad docking
newpy.write('\t\t\tcmd("mv %s Useless_dock"%(out_file))\n')
newpy.write('\t\t\tcmd("mv %s Useless_dock"%(arg[1]))\n')
newpy.write('\t\telif float(data[1]) >=-4.5:\n') # 0 > affinity >= -4.5 kcal/mol :: decent docking
newpy.write('\t\t\tcmd("mv %s Ok_dock"%(out_file))\n')
newpy.write('\t\t\tcmd("mv %s Ok_dock"%(arg[1]))\n')
newpy.write('\t\telse:\n') # -4.5 kcal/mol > affinity  :: good docking
newpy.write('\t\t\tcmd("mv %s Nice_dock"%(out_file))\n')
newpy.write('\t\t\tcmd("mv %s Nice_dock"%(arg[1]))\n')
newpy.write('\t\tcmd("rm %s"%(inp_file))\n')
newpy.close()

def multi_inp(cent_x = -18.287, cent_y = -9.034, cent_z = 7.361, max_dv = 5,
			  recpt = '6eqe.pdbqt', ligand = 'C8X.pdbqt', out_na = 'docked_C8X.pdbqt',
			  box_x = 15, box_y = 15, box_z = 15, exh = 10,
			  vina_home = '~/autodock_vina_1_1_2_linux_x86/bin'):

	coords_done = []
	sh = open('exec_vina.sh','w')
	sh.write('#!/bin/bash\n#\n#\n# Autor:Braga, B. C.\n# Email:bruno.braga@ufms.br\n#\n')
	sh.write('#'*30)
	sh.write('\nmkdir Useless_dock Ok_dock Nice_dock\n')

	i = 1
	while i <= multi:
		r_x = round(random.uniform(cent_x - max_dv, cent_x + max_dv), 3)
		r_y = round(random.uniform(cent_y - max_dv, cent_y + max_dv), 3)
		r_z = round(random.uniform(cent_z - max_dv, cent_z + max_dv), 3)
		r_coords = (r_x,r_y,r_z)
		if r_coords not in coords_done:
			coords_done.append(r_coords)
			f = open('config4vina_%d.inp'%i, 'w')
			f.write('receptor = %s\n'%recpt)
			f.write('ligand = %s\n\n'%ligand)
			f.write('out = R%d_%s\n\n'%(i,out_na))
			f.write('center_x = %f\n'%r_x)
			f.write('center_y = %f\n'%r_y)
			f.write('center_z = %f\n\n'%r_z)
			f.write('size_ = %.1f\n'%box_x)
			f.write('size_ = %.1f\n'%box_y)
			f.write('size_ = %.1f\n\n'%box_z)
			f.write('exhaustiveness = %d\n'%exh)
			f.close()

			sh.write('%s/vina --config config4vina_%d.inp --log Log_R%d.out\n'%(vina_home,i,i))
			sh.write('python3 verify_log.py Log_R%d.out\n'%i)
		else:
			continue
		i +=1
	sh.close()


if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv

	# Default keys
	exh       = 10
	box_x     = 15
	box_y     = 15
	box_z     = 15
	cent_x    = -18.287
	cent_y    = -9.034
	cent_z    = 7.361
	recpt     = '' #'6eqe.pdbqt'
	ligand    = '' #'C8X_BHET.pdbqt'
	out_na    = ''
	max_dv    = 5
	multi     = 10
	vina_home = '~/autodock_vina_1_1_2_linux_x86/bin'
	inst_only = False

	flags = ["&","-h","--help",'-vnh','-e','-dbox','-dcenter','-maxdev','-mult','-recpt','-lig','-out']

	# Flag verification
	for i in arg:
		if i[0] == '-' and i.lower() not in flags:
			print("Unkown Flag used: ", i)
			inst_only = True
			break

	i = 0
	while i < len(arg):
		if inst_only:
			break
		elif arg[i].lower() == "&":
			break
		elif arg[i].lower() == "-h" or arg[i].lower() == "--help":
			inst_only = True
			print("Copyright (C) 2021  Braga, B. C.\nThis program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions.\n")
			print("\nUsage:\n\t-h or --help\t\tPrints this message.\n")
			print("\t-vnh\t\tVina home directory.\n")
			print("\t-recpt\t\tReceptor file.\n\t-lig\t\tLigand file.\n")
			print("\t-out\t\tName of the pdbqt output file with all models vina propose. (Default: Docked_YOURLIGAND.pdbqt).\n")
			print("\t-dcenter\t\tXYZ-Coordinates of the center of the gridbox for docking.\n")
			print("\t-dbox\t\tDimensions XYZ of the gridbox.\n")
			print("\t-e\t\tExhaustiveness used. Default: 10.\n")
			print("\t-maxdev\t\tMaximum deviation from the gridbox center, for the random multiple random choices of gridboxes. Default: 5 angstrom.\n")
			print("\t-mult\t\tNumber of multiple docking tries. Default: 10.\n")
			print("\nExamples:\n\t$ python3 Multidock_Vna.py -recpt 6eqe.pdbqt -lig C8X_BHET.pdbqt\n")
			print("\n\t$ python3 Multidock_Vna.py -recpt 6eqe.pdbqt -lig C8X_BHET.pdbqt -out docked_C8X.pdbqt -vnh -dcenter -18.287 -9.034 7.361\n\t\tObs: '-vnh' used without any information if vina on your PC PATH\n")
			print("\n\t$ python3 Multidock_Vna.py -recpt 6eqe.pdbqt -lig C8X_BHET.pdbqt -out docked_C8X.pdbqt -dcenter -18.287 -9.034 7.361 -e 50 -maxdev 2.5 -mult 100\n")
			break

		# Key verifications
		elif arg[i].lower() == '-recpt':
			i+=1
			recpt = arg[i]
		elif arg[i].lower() == '-vnh':
			i+=1
			vina_home = arg[i]
		elif arg[i].lower() == '-lig':
			i+=1
			ligand = arg[i]
			out_na = 'Docked_%s'%ligand
		elif arg[i].lower() == '-out':
			i+=1
			out_na = arg[i]
		elif arg[i].lower() == '-dcenter':
			cent_x = float(arg[i+1])
			cent_y = float(arg[i+2])
			cent_z = float(arg[i+3])
			i+=3
		elif arg[i].lower() == '-maxdev':
			i+=1
			max_dv = float(arg[i])
		elif arg[i].lower() == '-mult':
			i+=1
			multi = int(arg[i])
		elif arg[i].lower() == '-dbox':
			box_x = float(arg[i+1])
			box_y = float(arg[i+2])
			box_z = float(arg[i+3])
			i+=3
		elif arg[i].lower() == '-e':
			i+=1
			exh = int(arg[i])
		i+=1

	if not inst_only:
		if '' not in [recpt, ligand]:
			multi_inp(cent_x=cent_x,cent_y=cent_y,cent_z=cent_z,max_dv=max_dv,recpt=recpt, ligand=ligand, out_na=out_na,box_x=box_x,box_y=box_y,box_z=box_z,exh=exh,vina_home=vina_home)
			cmd('sh exec_vina.sh')
		else:
			print('Receptor and/or ligand file not chosen!')
