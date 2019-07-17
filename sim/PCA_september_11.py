# -*- coding: utf-8 -*-
"""
Created on Fri May  5 09:48:51 2017

@author: amirbitran
"""

import numpy as np
import joblib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
#We will do PCA on coordinates data


#First, we need to obtain all coordinates, and align them relative to the reference

directories=['ADK_trunc26/replica_tight_spacing', 'ADK_trunc26/unfolding']   #directories containing results
temperatures=[['*.***'], ['*.***']]  #must have same number of sublists as directories. Each sublist corresponds to the temperatures we care about in that directory
modes=['replica', 'unfolding']  #what type of simulations are stored in each directory? We deal with them slightly differently for logistical reasons
#temperatures=[['*.***'], ['0.600']]
fileroot='adk_t26'
protein='ADK_trunc26'

combined_Data=[]  #will have data from all directories
combined_PDB_files=[]
ref_file=ref_file='{}/{}_0.900.0'.format(directories[0], fileroot)




for N, directory in enumerate(directories):
	temps=temperatures[N]
	if  os.path.exists("{}/Aligned_centered_snapshots.dat".format(directory)):
		[Data, mean, PDB_files]=joblib.load("{}/Aligned_centered_snapshots.dat".format(directory))   #data is centered
		Data=Data+mean  #we actually do not want the data to be centered, since we want to recenter at the end, once we have collected all of it, based on the combined mean of everything
		
		combined_Data.append(Data)
		combined_PDB_files.append([f for f in PDB_files])
		print('Loaded data from directory {}'.format(directory))
		
		#We need to double check that all the files we want are in Aligned_center_snapshots.dat

	else:  #To do: For unfolding simulations, chagne so that each temperature's centered files are saved separately
	    PDB_files=[]
	    for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) 
		
		#to do: SUPPRESS WARNINGS HERE!!!
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
	            
	    mean=np.mean(Data, axis=0)
	    #Data=Data-mean  #center the data
	    combined_Data.append(Data)
	    combined_PDB_files.append([f for f in PDB_files])
	    pickle.dump([Data-mean, mean, PDB_files], open("{}/Aligned_centered_snapshots.dat".format(directory),"wb"))   #We save the centered data, to be consistent with how things were done previously
	    print('Read and saved data from directory {}'.format(directory))   
	    
	    
		
		
		                                                                     


#Now we gotta do PCA
#First, we turn combined data into a numpy array

combined_Data=np.vstack(tuple(combined_Data))
combined_PDB_files=[f for sublist in combined_PDB_files for f in sublist]
print('Combined data incorporated into numpy array')

#Next, we center the combined data
combined_mean=np.mean(combined_Data, axis=0)
combined_Data=combined_Data-combined_mean
print('Combined data cetered')


pca=sklearn.decomposition.PCA(n_components=10)
#pca=sklearn.decomposition.PCA(n_components=np.shape(Data)[1])

pca.fit(combined_Data)  
projections=pca.fit_transform(combined_Data)
print('PCA performed on combined data')

pickle.dump([pca, projections, combined_mean, combined_PDB_files],open("{}/PCA_combined.dat".format(directories[0]), "wb"),protocol=pickle.HIGHEST_PROTOCOL )
print('Saved results. All done!')
#fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(projections[:,0], projections[:,1], projections[:,2])

