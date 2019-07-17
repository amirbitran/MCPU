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

directory='ADK/replica_tight_spacing3'   #CHANGE ON ODYSSEY
temps=['*.***']

if  os.path.exists("{}/Aligned_centered_snapshots.dat".format(directory)):
	[Data, mean, PDB_files]=joblib.load("{}/Aligned_centered_snapshots.dat".format(directory))
else:
    PDB_files=[]
    for temp in temps: PDB_files+=glob.glob('{}/adk_t26_{}*.*0'.format(directory, temp)) 
    
    ref_file='{}/adk_t26_0.900.0'.format(directory)
    reference=Bio.PDB.PDBParser().get_structure('ADK',ref_file).get_residues()
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
            s=Bio.PDB.PDBParser().get_structure('ADK',file)
            current=s.get_residues()
            structure=[]
            for res in current: structure.append(res['CA'])
            SI=Bio.PDB.Superimposer
            SI.set_atoms(SI,ref_structure, structure)
            SI.apply(SI,structure)
            coords=np.array([list(atom.get_coord()) for atom in structure])
            Data[n,:]=np.reshape(coords, [1, np.shape(coords)[0]*np.shape(coords)[1]])
    mean=np.mean(Data, axis=0)
    Data=Data-mean  #center the data

#Now we gotta do PCA




pca=sklearn.decomposition.PCA(n_components=10)
#pca=sklearn.decomposition.PCA(n_components=np.shape(Data)[1])

pca.fit(Data)  
projections=pca.fit_transform(Data)


pickle.dump([Data, mean, PDB_files], open("{}/Aligned_centered_snapshots.dat".format(directory),"wb"))                                                                            
pickle.dump([pca, projections, mean, PDB_files],open("{}/PCA_complete.dat".format(directory), "wb"),protocol=pickle.HIGHEST_PROTOCOL )

#fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(projections[:,0], projections[:,1], projections[:,2])

