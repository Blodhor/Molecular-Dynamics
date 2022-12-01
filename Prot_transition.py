def input_resname(arq="system.cpin"):
	f = open(arq,'r')
	temp = f.readlines()
	f.close()
	names = ''
	li = False
	for i in temp:
		if "RESNAME" in i:
			li = True
			for j in range(8,len(i)):
				if "'" not in i[j]:
					names += i[j]
			#print("0: ",i)
			continue
		if li:
			if "Residue:" in i:
				for j in range(1,len(i)):
					if "'" not in i[j]:
						names += i[j]
				#print("1: ",i)
			else:
				break
	data = names.split('\n')
	names = ''
	for i in data:
		names += i
	data = names.split(',')
	if "System:" in data[0]:
		res_id = data[1:]
	else:
		res_id = data
	if res_id[-1]=='':
		res_id.pop()
	return res_id

def cpout_read(arq='system.cpout', residues=[],total_time=200.0,compact=1):
	#total_time in nanosec
	#compact: 1==all data; 2==half; 3==a third...
	arq_name= arq.split('.')[0]
	f = open(arq,'r')
	temp = f.readlines()
	f.close()
	
	freq   = {}
	res_rt = {}
	for i in residues:
		freq[i]={}
		freq[i]['max']=0
		res_rt[i]=[] #state list
	
	count = 0
	for i in temp:
		if "Residue" in i:
			if count == len(residues):
				count = 0
			res_rt[residues[count]].append(i[-2])
			freq[residues[count]]['max'] +=1
			if "State: "+i[-2] not in freq[residues[count]]:
				freq[residues[count]]["State: "+i[-2]] = 1
			else:
				freq[residues[count]]["State: "+i[-2]] += 1
			count+=1

	for i in res_rt:
		tt = i.split(":") # Residue: TYR 59
		nm = tt[1].split()
		f = open("%s-%s_SvT.dat"%(arq_name,nm[0]+nm[1]),"w")
		f.write('Prot_State\tTime (ns)\n')
		saves = len(res_rt[i])
		step  = round(total_time/saves,4)
		time  = 0
		counter = 0 
		for j in res_rt[i]:
			time    += step
			counter += 1
			if counter%compact ==0:
				f.write('\t%s\t%.4f\n'%(j,time))
		f.close()
	
	freq_perc = {}
	for i in freq:
		da = i.split(":")
		freq_perc[da[1][1:]] = {}
		for j in freq[i]:
			if j !="max":
				freq_perc[da[1][1:]][j] ="%.2f%%"%round(freq[i][j]*100.0/freq[i]['max'],2)
	
	f = open("%s-Freq.README"%arq_name,"w")
	for i in freq_perc:
		f.write('%s: '%i)
		for j in freq_perc[i]:
			f.write("("+j+",%s),"%freq_perc[i][j])
		f.write("\n")
	f.close()
	return freq_perc

if __name__ == "__main__":
	from os import system as cmd
	import sys
	arg = sys.argv

	#default lines
	inst_only  = False
	cpin_name  = ''
	cpout_name = ''
	compact    = 500
	dur        = 200.0 #ns
	flags = ["&","-h","--help",'-cpin','-cpout','-compact','-duration']

	#pode muda a vontade 'path','cpin','cpout' e 'Default'
	#  que nao altera nada se tiver flag de entrada
	Default = True
	path     = 'Dev_cpoutAnalyser/'
	enzyme   = ['Nat','D206EH237K'][0]
	md_id    = {1:1,3:3}[3]
	inp_name = '%s_CpH7MD%d'%(enzyme,md_id)
	cpin     = '%s.cpin'%inp_name
	cpout    = '%s.cpout'%inp_name
	
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
			print("\t-compact\t\tCompression factor.Ex: 1 means all data, 10 means a tenth.\n")
			print("\t-duration\t\tDuration of your production stage in nano seconds.\n")
			print("\t-cpin\t\tYour System cpin file.\n")
			print("\t-cpout\t\tYour System production stage cpout file.\n")
			print("\t-Ex:\n\t\tpython.exe Prot_transition.py -cpin system.cpin -cpout Prod.cpout -duration 200 -compact 100\n")
			break

		# Key verifications
		elif arg[i].lower() == '-compact':
			i+=1
			compact = int(arg[i])
		elif arg[i].lower() == '-duration':
			i+=1
			dur = float(arg[i])
		elif arg[i].lower() == '-cpin':
			Default=False
			i+=1
			cpin_name = arg[i]
		elif arg[i].lower() == '-cpout':
			Default=False
			i+=1
			cpout_name = arg[i]
		i+=1

	if Default:
		cpin_name  = path+cpin
		cpout_name = path+cpout

	if not inst_only:
		res = input_resname(cpin_name)
		cpout_read(arq=cpout_name,residues=res,total_time=dur,compact=compact)
