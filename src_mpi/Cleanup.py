import os
import numpy as np

files=os.listdir()


for f in files:
	fields=f.split('.')
	if len(fields)==3 and fields[3]!='log':
		number=int(fields[2])
		if np.mod(number, 500000)!=0:
			os.remove(f)
	