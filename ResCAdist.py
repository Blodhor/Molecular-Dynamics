# Distance betwwen two atoms
def dist(atm1=(1,-2,3), atm2=(0,5,-6)):
	x1, y1, z1 = atm1
	x2, y2, z2 = atm2
	r = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
	return round(r,2)

def read_ca_coords(res=['178','209'], pdb='bla.pdb'):
	f = open(pdb)
	temp = f.readlines()
	f.close()

	ca = []
	count = 0
	for i in temp:
		if count == len(res):
			break
		if 'ATOM ' in i:
			data = i.split()
			if data[2] == 'CA' and data[4] in res:
				count += 1
				print(i)
				ca.append( (float(data[5]),float(data[6]),float(data[7])) ) 
	return ca

def Main(arg=['ResCAdist.py','file.pdb','178','209']):
	atms = read_ca_coords(res=[arg[-2],arg[-1]],pdb=arg[1])
	print("Distance =", dist(atm1=atms[0],atm2=atms[1]))

def Help(arg=['ResCAdist.py', '-h']):
	print("Please run this code as follows:\n")
	print("\t$ pythonCompilerV3+ %s PDBFILE.pdb res1 res2"%arg[0])

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv
	
	try:
		if len(arg) == 4 and '.pdb' in arg[1]:
			Main(arg=arg)
		else:
			Help()
	except:
		Help()