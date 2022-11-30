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
			if "RESSTATE" not in i:
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

def cpout_read(arq='system.cpout', residues=[],total_time=200):
	#total_time in nanosec
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
		step  = round(1.0*total_time/saves,4)
		time  = 0 
		for j in res_rt[i]:
			time += step
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

path = 'Dev_cpoutAnalyser/'
cpin  = 'D206EH237K_CpH7MD1.cpin'
cpout = 'D206EH237K_CpH7MD1.cpout'
res = input_resname(path+cpin)
print(cpout_read(path+cpout,res))
