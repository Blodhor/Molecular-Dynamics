def Create_averaged_data(files=[], new_file='RMSF_mean.dat'):
	X=[]
	Y=[]
	intx = False
	if 'rmsd' in files[0]:
		new_file='RMSD_mean.dat'
		intx = True
	elif 'radgyr' in files[0]:
		new_file='RADGYR_mean.dat'
		intx = True
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
			#it is recommended to use files with the same size,
			#  but it will still work with the script below
			#  (if their size difference is not too big)
			t = len(Y)-len(y)
			if t>0:
				#y<Y
				cp = y[len(y)-t:]
				for j in cp:
					x.append(x[-1]+1)
					y.append(j)
				print("Fixed the size difference of ", len(cp))
			elif t<0:
				#Y<y
				cp = Y[len(Y)+t:]
				for j in cp:
					X.append(X[-1]+1)
					Y.append(j)
				print("Fixed the size difference of ", len(cp))
					
			for j in range(len(y)):
				Y[j] += y[j]

	for i in range(len(X)):
		if intx:
			sx= "\t%d"%X[i]
		else:
			sx= "%.3f\t"%X[i]
		new.write("%s\t%.4f\n"%(sx, Y[i]))
	new.close()
	
	return 0


if __name__ == "__main__":
	import sys
	arg = sys.argv

	#print(arg)
	if len(arg) > 2:
		Create_averaged_data(files=arg[1:])
	else:
		print("You must give TWO or more RMSF files on the 2 colunm style of .dat!")
