def Create_averaged_data(files=[], new_file='RMSF_mean.dat'):
	X=[]
	Y=[]
	new = open(new_file,'w')
	n = len(files)
	if n<2:
		print("I need more files!")
		return -1
	count =0
	for i in files:
		f = open(i,'r')
		temp = f.readlines()
		f.close()
		x=[]
		y=[]
		for j in temp:
			if '#' in j:
				if count==0:
					count +=1
					new.write(j)
				continue
			data = j.split()
			if len(data) == 2:
				x.append( float(data[0]) )
				y.append( float(data[1])/n )
		
		if len(Y)==0:
			for j in range(len(y)):
				X.append(x[j])
				Y.append(y[j])
		else:
			for j in range(len(y)):
				Y[j] += y[j]

	for i in range(len(X)):
		new.write("%.3f\t\t%.4f\n"%(X[i], Y[i]))
	new.close()
	
	return 0


if __name__ == "__main__":
	import sys
	arg = sys.argv

	#print(arg)
	if len(arg) > 2:
		Create_averaged_data(files=arg[1:], new_file='RMSF_mean.dat')
	else:
		print("You must give TWO or more RMSF files on the 2 colunm style of .dat!")
