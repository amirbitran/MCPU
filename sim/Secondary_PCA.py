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
#We will do PCA on coordinates data for those files whose projection into the first
#principal component in the original PCA is negative


#First, we need to obtain all coordinates, and align them relative to the reference

directory='ADK/unfolding'   #CHANGE ON ODYSSEY
temp='0.900'
[pca, projections, mean, unfolding_PDB_files]=joblib.load('{}/PCA.dat'.format(directory))
indices=np.where(projections[:,0]<0)[0]  #which indices correspond to files that have a negative first principal component projection?

PDB_files=[unfolding_PDB_files[i] for i in indices]

ref_file='{}/adk_0.900_21.0'.format(directory)
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




#Now we gotta do PCA


#First, we center the data
mean=np.mean(Data, axis=0)

centered_data=Data-mean


pca2=sklearn.decomposition.PCA(n_components=3)

pca2.fit(Data)  
projections=pca2.fit_transform(Data)


pickle.dump([pca2, projections, mean, PDB_files],open("{}/Secondary_PCA.dat".format(directory), "wb"),protocol=pickle.HIGHEST_PROTOCOL )

#fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(projections[:,0], projections[:,1], projections[:,2])

