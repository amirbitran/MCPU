import numpy as np
import os
import glob

directory='ISPA_trunc87/unfolding'
fileroot='ispa_t87'

def Get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times
    
    

print('Obtaining list of PDB files')


temp='*.***'

#All the PDB files we currently have...gross..wayy too many
PDB_files=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) #USE ON  Odyssey


print('{} PDB files at the moment...gonna take a while'.format(len(PDB_files)))

for f, file in enumerate(PDB_files):
	if np.mod(f,500)==0:
		print(f)
	time=Get_times([file])[0]
	if np.mod(time, 500000)!=0:
		os.remove(file)