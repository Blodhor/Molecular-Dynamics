def reformat(name='nome',new='novo-nome',shift=28):
	default_counter = {'HIP': 'HIS', 'AS4': 'ASP', 'GL4': 'GLU'}
	f = open(name,'r')
	temp = f.readlines()
	f.close()
	
	for i in range(len(temp)):
		if 'ATOM' in temp[i]:
			z = temp[i].split()
			for j in default_counter:
				if j in z[3]:
					z[3]= default_counter[j]
					break
			z[4] = str(int(z[4])+shift)
			temp[i] = z[0]+'    '+z[1]
			if len(z[2]) == 1:
				temp[i] += '    '+z[2]+'    '+z[3]
			elif len(z[2]) == 2:
				temp[i] += '   '+z[2]+'   '+z[3]
			elif len(z[2]) == 3:
				temp[i] += '  '+z[2]+'  '+z[3]
			else:
				temp[i] += ' '+z[2]+' '+z[3]
			for j in z[4:]:
				temp[i] += '    '+j
			temp[i]+='\n'

	f = open(new,'w')
	for i in temp:
		f.write(i)
	f.close()

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv

	#default lines
	inst_only  = False
	flags = ["&","-h","--help",'-i','-out','-rs']

	name='nome'
	new='new.pdb'
	shift = 0#28
	
	# Flag verification
	for i in arg:
		if i[0] == '-':
			try:
				if i.lower() not in flags:
					print("Unknown Flag used: ", i)
					inst_only = True
					break
			except:
				print("Unknown input comand, please check the --help option!")

	i = 0
	while i < len(arg):
		if inst_only:
			break
		elif arg[i].lower() == "&":
			break
		elif arg[i].lower() == "-h" or arg[i].lower() == "--help":
			inst_only = True
			print("\nUsage:\n\t-h or --help\t\tPrints this message.\n")
			print("\t-i\t\tPDB file (AMBER format).\n")
			print("\t-out\t\tName for the normal formatted PDB file.\n")
			print("\t-rs\t\tShift on the residue id (Default: 0).\n")
			print("\t-Ex:\n\t\tpython.exe Reformat_amberPDB.py -i MD.pdb -out MD_reform.pdb -rs 28\n")
			break

		# Key verifications
		elif arg[i].lower() == '-rs':
			i+=1
			shift = int(arg[i])
		elif arg[i].lower() == '-i':
			i+=1
			name = arg[i]
		elif arg[i].lower() == '-out':
			i+=1
			new= arg[i]
		i+=1

	if not inst_only:
		reformat(name=name,new=new,shift=shift)