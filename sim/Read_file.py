                                                                             
import numpy as np
#import pymbar # for MBAR analysis
#from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
import glob


##############################################################################
# First, read the file
##############################################################################

directory='CMK_trunc98/MultiUmbrella5'
temperature='All'
variables=['natives', 'energy', 'rmsd']
#variables=['potnl', 'sctor', 'hbond', 'Aro', 'torsion']
#variables=['energy']
step_multiples_to_read=1 #only read steps if they are a multiple of this number
min_step = 0 #only read steps that are greater than this number (inclusive)
max_step=np.inf  #Only read steps that are less than this number

filename='All_Data.dat'


def get_temp(filename):
	splitline=filename.split(sep='/')
	split2=splitline[2].split('_')
	#print(split2)
	while split2[1][0] not in '0123456789':
		del split2[1]
	temp=float(split2[1][0:5])
	print(temp)
	return temp


def read_file(PDB_files, variables):
	"""
	PDB_files actually means log_files...lol sorry
	
	variables is a list of variables that you want to read from the log files, but they need to be called
	exactly what they are called in the first line of the log file
	
	for instance, 'energy', 'natives', etc...
	
	Returns a 3D array data where data[i,j,k] corresponds to log_file i, time j within that log file, and variable k (in case you care about multiple variables like energies, natives, etc)
	
	"""
	data=[]
	lens=[]
	variable_indices=[]
	times=[]
	temperatures=[]
	setpoints=[]
	for filecounter, filename in enumerate(PDB_files):
		step_index=0
		print("Reading file {}".format(filename))    
		openfile=open(filename)
		data.append([])
		#energies.append([])
		#contacts.append([])
		#rmsd.append([])
		 #temperatures.append(float(temp))
		#       file 1   variable 1  times     variable 2  times
		#data=[ [        [  x1, x2, x3  ... ],  [y1, y2, y3...   ],...   ]   ]
		for line in openfile.readlines():
			line=line.rstrip('\n')
			if len(line)>0:
				entries=line.split()
				if 'step #' in line:
					fields = ['step'] + line.split()[2:]
					temperature_index=fields.index('temp')
					if 'setpoint' in fields:
						setpoint_index=fields.index('setpoint')
					else:
						setpoint_index=np.nan
					for variable in variables:
						variable_indices.append(fields.index(variable))
						data[filecounter].append([])
					#print(variable_indices)
				if entries[0]=='STEP':
					if np.mod(int(entries[1]), step_multiples_to_read)==0 and int(entries[1])>=min_step and int(entries[1])<max_step:
						step_index+=1
						if filecounter==0:
							times.append(int(entries[1]))
							#print(entries[variable_indices[1]+1])
						if step_index==1:  #learn what reporter values we currently have...only need to do this once per log file
							temperatures.append(float(entries[temperature_index+1][0:5]))
							if 'setpoint' in fields:
								setpoints.append(float(entries[setpoint_index+1]))
							else:
								setpoints.append(0)			
						for v, variable in enumerate(variables):
							data[filecounter][v].append(float(entries[variable_indices[v]+1]))
		

		lens.append(len(data[filecounter][0]))  
		data[filecounter]=np.array(data[filecounter])

		#if filecounter==0:
		#	print(data[0][1,:])
		x=np.zeros((1, len(data[filecounter][0]), len(data[filecounter])))
		for v in range(len(variables)):
			x[0,:,v]=data[filecounter][v,:]
		data[filecounter]=x
		#if filecounter==0:
		#	print(data[filecounter][0,:,1])		
	
	nonzero_lengths = [i for i in range(len(lens)) if lens[i]>0]
	
	data = [x for i, x in enumerate(data) if i in nonzero_lengths]
	lens = [l for i, l in enumerate(lens) if i in nonzero_lengths]
	
	data=np.vstack((x[:, 0:min(lens), :] for x in data))
	#print(data[0,:,1])
	#We have created an array data that is files (i.e. conditions) by timepoints by variables
	return data, temperatures, setpoints, np.array(times) 
    





All_files=glob.glob('{}/*.log'.format(directory))



log_files=[]

if temperature!='All':
	for file in All_files:
		#print(file)
		if get_temp(file)==temperature:
			print(file)
			log_files.append(file)
else:
	log_files=All_files
	
		#print('{}: {}'.format(file, get_temp(file)))
		#if get_temp(file) not in temperatures:
			#print('{} should not be there'.format(file))
			#PDB_files.remove(file)
			
			
#print (PDB_files)  
print('Reading files')        
data, temperatures, setpoints, times=read_file(log_files, variables)
#print(data[0,:,1])


print("Saving data")

pickle.dump([data, variables, log_files, temperatures, setpoints, times],open("{}/{}".format(directory, filename), "wb"), -1)
