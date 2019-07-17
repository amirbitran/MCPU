# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:27:32 2017
@author: amirbitran

Computes a mean contact map for every cluster
All you need to input are a list of PDB files and a list of labels which tell you the cluster assignment for each PDB file
"""

import numpy as np
import joblib
import matplotlib.pyplot as plt
import scipy
from scipy import linalg
import Analyze_structures
import pickle
import fnmatch



directory='ADK/replica_tight_spacing3'
atom='CA'  #we care about the alpha carbons
traj='***_*.***.'  #We care about all temperatures but specifically, replica exchange simulations
save_suffix=None
#save_suffix='lowtemps'  #how we will call the results
traj='0.6**.'

[gmix,labels, PDB_files]=joblib.load('{}/Gaussian_labels.dat'.format(directory))
        
def lookup_files(label,labels, PDB_files, traj='All'):
    """
    returns all files that are assigned to a given label   
    You can use traj to specify that you only care about files in a given trajectory
    The defalt is 'All' (all trajectories are considered )
    """
    files=[f for n,f in enumerate(PDB_files) if labels[n]==label ]
    if traj!='All':
         files = [f for f in files if fnmatch.fnmatch( f, '*_{}*'.format(traj))]
    return files
         
            
reference_coords, resis=Analyze_structures.read_PDB(PDB_files[0], atom)
reference_matrix=Analyze_structures.compute_contacts_matrix(reference_coords)


#we will now make a matrix where the ith page corresponds to the mean contacts map of the ith cluster
cluster_contact_maps=np.zeros((len(reference_matrix), len(reference_matrix), len(set(labels))))


for n in set(labels):
	print('Reading files in cluster{}'.format(n))
	if traj!='All':
		cluster=lookup_files(n, labels, PDB_files,traj)   #all files corresponding to this label
	else:
		cluster=lookup_files(n, labels, PDB_files, '*')
	local_maps=np.zeros((len(reference_matrix), len(reference_matrix), len(cluster)))   #set of contact maps for this particular cluster
	c_term_contact=0
	for i, file in enumerate(cluster):
		coords, resis=Analyze_structures.read_PDB(file, atom)
		local_maps[:,:,i]=Analyze_structures.compute_contacts_matrix(coords, thresh=8, min_seq_separation=4)  
		#if np.sum(local_maps[200:215,168:183,i])>1: c_term_contact+=1	#check whehter there is a C-terminal kinetic trap--purely for informative purposes
	cluster_contact_maps[:,:,n]=np.mean(local_maps, axis=2)
	#print('The fraction of snapshots assigned to label {} containing a C terminal kinetic trap is {}'.format(n, c_term_contact/len(cluster)))

if save_suffix==None:	
	pickle.dump([PDB_files, labels, cluster_contact_maps],open("{}/ClusterContactMaps_{}.dat".format(directory, atom), "wb"))
else:
	pickle.dump([PDB_files, labels, cluster_contact_maps],open("{}/ClusterContactMaps_{}_{}.dat".format(directory, atom, save_suffix), "wb"))
