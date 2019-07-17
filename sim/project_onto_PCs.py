# -*- coding: utf-8 -*-
"""
Created on Mon May  8 10:05:08 2017
This script projects snapshots from a given set of trajectories onto principal components
obtained by doing PCA previously on some other set of trajectories
@author: amirbitran
"""


import numpy as np
import joblib
from scipy import linalg
import sklearn
import glob
import Bio
from Bio import PDB
import sklearn
from sklearn import decomposition
import pickle
import warnings






#First, we need to obtain all coordinates, and align them relative to the reference

directory='ALP/MultiUmbrella_fixed2'    #contains snapshots you want to project  
PCA_directory='ALP/MultiUmbrella_orig2'   #contains existing PCA results
root = 'alp'

#first, we list all the files.

PDB_files=[glob.glob('{}/{}_*.*0'.format(directory, root))] 


#we make a list of lists, where each sublist contains all the files at a given temperature

mintemp=0.600
maxtemp=0.600
spacing=1.00
temps=np.arange(mintemp, maxtemp+0.00001, spacing)   #add 0.0001 so that the upper bound (maxtemp) is included 


#max_mcstep=199000000
#PDB_files=[]
#for t in temps:
#    t = str(t)
#    while len(t)<5: t+='0'  #all temps should be of the format 0.450, with three digits after decimal
    #PDB_files.append([ '{}/adk_{}.{}'.format(directory, t, int(n)) for n in np.arange(0, max_mcstep, 1) ])
#    PDB_files.append( glob.glob('{}/adk_{}_*.*0'.format(directory, t)) )    

#Now, we load the reference file--every file will be aligned to this one

#ref_file='{}/{}_0.600.0'.format(PCA_directory, root)
ref_file = 'ALP/MultiUmbrella_orig/alp_0.600_Emin.pdb'


reference=Bio.PDB.PDBParser().get_structure('Structure',ref_file).get_residues()
ref_structure=[]
for res in reference: ref_structure.append(res['CA'])
ref_coords=np.array([list(atom.get_coord()) for atom in ref_structure])
    

    

#load PCA components from unfolding simulations
[pca, projections, mean_source_snap, source_PDB_files]=joblib.load('{}/PCA_combined.dat'.format(PCA_directory))



#We create the array which will contain all the data
#Each page will hold one temperature
#Within that page, rows will be files at that temperature, and columns will be dimensions (of whcih there are n_atoms*3, for x,y,z)
#Data=np.zeros((max([len(PDB_files[i])for i in range(len(PDB_files))]),np.shape(ref_coords)[0]*np.shape(ref_coords)[1], len(temps) ))  
#not all temperatures will have the same number of files, so we preallocate the Data array to have the size for the largest temp
#There will be some empty 0 entires as a result

#now we obtain the components themselves
components=pca.components_
n_comps=np.shape(components)[0]


#file=PDB_files[1]
for t in range(len(temps)):
    print('Working on temperature {}'.format(temps[t]))
    Data=np.zeros((  len(PDB_files[t]) ,  np.shape(ref_coords)[0]*np.shape(ref_coords)[1]   ))
    for n, file in enumerate(PDB_files[t]):
        with warnings.catch_warnings():  #suppress all warnings
            warnings.simplefilter("ignore")
            #print('Analyzing {}'.format(file))
            s=Bio.PDB.PDBParser().get_structure('Structure',file)
            current=s.get_residues()
            structure=[]
            for res in current: structure.append(res['CA'])
            #orig_coords=[list(atom.get_coord()) for atom in structure]
            
            SI=Bio.PDB.Superimposer
            SI.set_atoms(SI,ref_structure, structure)
            SI.apply(SI,structure)
            coords=np.array([list(atom.get_coord()) for atom in structure])
            Data[n,:]=np.reshape(coords, [1, np.shape(coords)[0]*np.shape(coords)[1]])
    centered_data=Data-mean_source_snap  #center the data
    
    #project data onto components
    replica_projections=np.dot(centered_data, np.transpose(components))
    
    
#(n_files x n_dimensions) * (n_dimensions x n_components) = (n_files x n_components)

######f1#####   #  #  #  
######f2#####   #  #  #  
######f3#####   v1 v2 v3 
######f4#####   #  #  #
                #  #  #
    temp = str(temps[t])
    while len(temp)<5: temp+='0'
    pickle.dump([pca, replica_projections, mean_source_snap, PDB_files[t]],open("{}/PCA_combined.dat".format(directory, temp), "wb"),protocol=pickle.HIGHEST_PROTOCOL )

    
           
    



