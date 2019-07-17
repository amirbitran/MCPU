# -*- coding: utf-8 -*-
"""
Created on Fri May  5 09:48:51 2017
Obtains coordinates for all data, and combines data into a single matrix, and does PCA

The code is written in such a way that you can choose whcih temperatures you want data from
But nonetheless, you keep a repository of all the data with which you've EVER worked in Aligned_centered_snapshots.dat


@author: amirbitran
"""

import numpy as np
import joblib
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy import linalg
import sklearn
import glob
import Bio
from Bio import PDB
import sklearn
from sklearn import decomposition
import pickle
import warnings
import os
import sys
#We will do PCA on coordinates data



#Send print output to a special file
#:q
#sys.stdout = open('PCA_t26_out.txt', 'w')

#First, we need to obtain all coordinates, and align them relative to the reference

#directories=['ALP/MultiUmbrella_orig', 'ALP/MultiUmbrella_fixed', 'ALP/MultiUmbrella_nolocal']   #directories containing results
directories=['ALP/MultiUmbrella_orig2', 'ALP/MultiUmbrella_fixed2', 'ALP/MultiUmbrella_nolocal2']
#directories=['MAP/replica_tight_spacing']
#temperatures=[['*.***']]  #must have same number of sublists as directories. Each sublist corresponds to the temperatures we care about in that directory
temperatures=[['*.***'], ['*.***'], ['*.***']]
fileroot='alp'
protein='ALP'
n_components=10

combined_Data=[]  #will have data from all directories
combined_PDB_files=[]
#ref_file=ref_file='{}/{}_0.600.0'.format(directories[0], fileroot)
ref_file = 'ALP/MultiUmbrella_orig/alp_0.600_Emin.pdb'

#which residues do we care about in our PCA?

residues='All'
#residues=np.arange(70,264,1)  
def Get_aligned_files(PDB_files, protein, ref_file):
	
	#to do: SUPPRESS WARNINGS HERE!!!
	with warnings.catch_warnings(): #suppress all warnings
		warnings.simplefilter("ignore")
		reference=Bio.PDB.PDBParser().get_structure(protein,ref_file).get_residues()
		ref_structure=[]
		for res in reference: ref_structure.append(res['CA'])
		ref_coords=np.array([list(atom.get_coord()) for atom in ref_structure])
    
    
	#Data will be our matrix of observations by features
	#Each coordinate of each atom will be a feature, so the total nubmer of features
	#is 3*n_atoms
	Data=np.zeros((len(PDB_files),np.shape(ref_coords)[0]*np.shape(ref_coords)[1] ))  
    

    

	#file=PDB_files[1]
	for n, file in enumerate(PDB_files):
	    with warnings.catch_warnings():  #suppress all warnings
	        warnings.simplefilter("ignore")
	        s=Bio.PDB.PDBParser().get_structure(protein,file)
	        current=s.get_residues()
	        structure=[]
	        for res in current: structure.append(res['CA'])
	        SI=Bio.PDB.Superimposer
	        SI.set_atoms(SI,ref_structure, structure)
	        SI.apply(SI,structure)
	        coords=np.array([list(atom.get_coord()) for atom in structure])
	        Data[n,:]=np.reshape(coords, [1, np.shape(coords)[0]*np.shape(coords)[1]])
	return Data


for N, directory in enumerate(directories):
	temps=temperatures[N]
	if  os.path.exists("{}/Aligned_centered_snapshots.dat".format(directory)):
		[Data, mean, PDB_files]=joblib.load("{}/Aligned_centered_snapshots.dat".format(directory))   #data is centered
		Data=Data+mean  #we actually do not want the data to be centered, since we want to recenter at the end, once we have collected all of it, based on the combined mean of everything
		print('Loaded data from directory {}'.format(directory))
		
		original_Data=Data  #copy this variable since we will need it later in its intact form
		original_files=PDB_files
		
		#First thing: We see if there are any files and datapoints in here that we don't need
		#needed_files are all the files we need
		needed_files=[]
		for temp in temps: needed_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) 
				
		indices_to_remove=[i for i in range(len(PDB_files)) if PDB_files[i] not in needed_files]
		
		
		Data=np.delete(Data, indices_to_remove, 0) #removes the rows (hence along axis "0") indicated by the above indices
		PDB_files = [PDB_files[i] for i in range(len(PDB_files)) if i not in indices_to_remove]
		
		
		combined_Data.append(Data)
		combined_PDB_files+=[f for f in PDB_files]
		
		
		#We need to double check that all the files we want are in Aligned_center_snapshots.dat
		#If some are missing, we should add them using Get_aligned_files function
		
		
		#Now check if all files in PDB files are in needed_files
		#missing_files is a list of files that are in needed_files but not in PDB_files
		
		missing_files=[f for f in needed_files if f not in PDB_files]
		
		#Feed missing_files into Get_aligned_files to obtain the missing coordinates, and append accordingly
		
		if len(missing_files)>0:
			print('Some files are missing!')
			missing_data=Get_aligned_files(missing_files, protein, ref_file)
			print('Obtained data for additional PDB files')
			combined_Data.append(missing_data)
			combined_PDB_files+=[f for f in missing_files]
			
			Data_to_save=np.vstack((original_Data, missing_data))  #we combine the original data (without having removed any files) with the data we accumulated during this run, and save everything, so that our repository grows
			new_mean=np.mean(Data_to_save, axis=0)
			files_to_save=original_files + missing_files
			pickle.dump([Data_to_save-new_mean, new_mean, files_to_save], open("{}/Aligned_centered_snapshots.dat".format(directory),"wb"))
			
		
			
		
		
		
		
		
		

	else:  #There is no saved centered data in this directory
		PDB_files=[]
		for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) 
		
	
		Data= Get_aligned_files(PDB_files, protein, ref_file)
		mean=np.mean(Data, axis=0)
		combined_Data.append(Data)
		combined_PDB_files+=[f for f in PDB_files]
		pickle.dump([Data-mean, mean, PDB_files], open("{}/Aligned_centered_snapshots.dat".format(directory),"wb"))   #We save the centered data, to be consistent with how things were done previously
		print('Read and saved data from directory {}'.format(directory))   
	    
	    
		
		
		                                                                     


#Now we gotta do PCA
#First, we turn combined data into a numpy array

combined_Data=np.vstack(tuple(combined_Data))
#combined_PDB_files=[f for sublist in combined_PDB_files for f in sublist]
print('Combined data incorporated into numpy array')

#Next, we center the combined data
combined_mean=np.mean(combined_Data, axis=0)
combined_Data=combined_Data-combined_mean
print('Combined data cetered')





pca=sklearn.decomposition.PCA(n_components=n_components)
#pca=sklearn.decomposition.PCA(n_components=np.shape(Data)[1])



#we figure out what indicies we want to do PCA on, depending on which residues we care about
if residues=='All':
	pca.fit(combined_Data)
	projections=pca.fit_transform(combined_Data)
else: 
	indices=np.arange(3*(residues[0]-1), 3*residues[-1],1)  
	pca.fit(combined_Data[:,indices])
	projections=pca.fit_transform(combined_Data[:, indices])
print('PCA performed on combined data')

if residues=='All':
	pickle.dump([pca, projections, combined_mean, combined_PDB_files],open("{}/PCA_combined.dat".format(directories[0]), "wb"),protocol=pickle.HIGHEST_PROTOCOL )
else: 
	pickle.dump([pca, projections, combined_mean, combined_PDB_files],open("{}/PCA_{}_through_{}.dat".format(directories[0], residues[0], residues[-1]), "wb"),protocol=pickle.HIGHEST_PROTOCOL )
print('Saved results. All done!')
#fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(projections[:,0], projections[:,1], projections[:,2])

