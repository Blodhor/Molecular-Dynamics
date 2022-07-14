#extraindo dados do mmpbsa_decomp.dat
def transfo(old_data = "m.txt",
new_data = "pH7-MD1_PB-BindEnerg.dat",
linhas=range(7,273)):
	f = open(old_data,'r')
	m = f.readlines()
	f.close()

	data =[]
	lines = [i for i in linhas] #total energy decomp para IsPETase+BHET

	for i in lines:
		ll = m[i].split(',')
		# residue + lig(L)/recp(R) , Total binding energy +- std. deviation
		data.append( (ll[0]+"\t"+ll[1][0], ll[-3]+" +- "+ll[-2]) ) 

	f = open(new_data,"w")
	f.write("Total = internal + van der waals + eletrostatic + polar solvation + non-polar solvation\n\n")
	f.write("Res\tLoc\tTotal +- Std. Deviation\n")
	for i in data:
		f.write(i[0]+"\t"+i[1]+"\n")

	f.close()

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv

	flags = ["&","-h","--help",'-i','-out']

	# Flag verification
	for i in arg:
		if i[0] == '-':
			try:
				if type(float(i)) == type(2.3):
					#pra aceitar numero sozinho
					continue
			except ValueError:
				if i.lower() not in flags:
					print("Unkown Flag used: ", i)
					inst_only = True
					break
	#default lines
	lini = 7
	lfinal = 273


	i = 0
	while i < len(arg):
		if inst_only:
			break
		elif arg[i].lower() == "&":
			break
		elif arg[i].lower() == "-h" or arg[i].lower() == "--help":
			inst_only = True
			print("\nUsage:\n\t-h or --help\t\tPrints this message.\n")
			print("\t-i\t\tFinal Decomp data for PB in MMPBSA.py.\n")
			print("\t-out\t\tName of the output file with only the residue and total energy.\n")
			print("\t-lines\t\tInput file lines for residues total decomp. Please indicate the numbers of the initial res. line and the empty line after all data separated by comma.\n")
			print("\t-Ex:\n\t\tpython.exe -i YourInput -out NameYourOutput -lines 8,274\n")
			print("\t\tPS: Please do not use spaces in your file's names!")
			break

		# Key verifications
		elif arg[i].lower() == '-i':
			i+=1
			input_name = arg[i]
		elif arg[i].lower() == '-out':
			i+=1
			output_name = arg[i]
		elif arg[i].lower() == '-lines':
			i+=1
			temp  = arg[i].split(',')
			lini  = int(temp[0])-1
			lfinal= int(temp[1])-1
		i+=1

	if not inst_only:
		transfo(old_data=input_name,new_data=output_name,linhas=range(lini,lfinal))