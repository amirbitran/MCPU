#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 08:13:37 2018

@author: amirbitran
"""


import numpy as np
import os
import joblib
import itertools
import joblib
import Analyze_structures
import LoopCluster
import glob
import copy as cp
import natsort
import fnmatch
import copy as cp
import ClusterSubstructures



protein='MARR_trunc44' #eax. MARR or MARR_trunc44
fileroot='marr_t44' #ex. marr or marr_t44
mode = 'binary'
ID_nonnatives = False  #Do we care to count how many nonnative substructures are present in each snapshot?


#If we want to specify the substructures in advance, we include the score file containing those substructures in the following path
#Note that this could be the substructures for a larger version of this protein (ex. for MARR_trunc44, you may specify the Substructure_file for full MARR)
#In that case, you only keep substructures/portions of substructures that involve amino acids present in this current (shorter) protein

#Alternatively, you are allowed to use substructures for a smaller protein (ex. full MARR) to analyze a larger protein (ex. MARR dimer)
#In that case, you only use those amino acids of the larger protein that are also present in the smaller full protein

Substructure_file='MARR/MultiUmbrella3/Substructure_scores.dat'
#Substructure_file='CMK/unfolding3/Substructure_scores.dat'
###at what distance do we define a contact? ###
d_cutoff=7.8   #6.5 for FABG and CMK, 7.8 for MARR  
trap_d_cutoff = 8.5  #Cutoff for identifying nonnative substructures,  more lenient than native substructures...8.5 for MARR, let's try 7.5 for FABG
##########

directories=['{}/MultiUmbrella3'.format(protein)]  #what directories we analyze
#native_file='MARR/files/marr_0.100_Emin.pdb'
native_file='{}/files/{}_0.100_Emin.pdb'.format(protein, fileroot)

temperatures=[['*.***']]
#temperatures=[[ '0.650', '0.675', '0.700']] #previously [['*.***']]


output_filename='Substructure_contacts.dat'

### The following parameters are relevant if you do not specify a substructure file ####

min_clustersize=7 #7 for MARR, 7 for FABG
contact_sep_thresh=5 #7 for FABG, 3 for MARR, 5 for CMK





def Score_snapshot(snapshot, substructures, mode, d_cutoff, atom='CA', min_seq_separation=8,  ):
    """
    Assigns a set of scores for a snapshot
    This can be one of two things:
    
    1. If mode = 'distances', 
    the ith score tells you what is the average distnace between pairs of residues residues that participate in the ith substructure, in this snapshto

	2. If mode = 'binary'
	the ith score tells you how many native contacts that belong to the ith substructure have formed
    """
    if mode not in ['distances', 'binary']:
    	print("Error! Value for mode must be either 'distances' or 'binary'")
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    distances=Analyze_structures.compute_contacts_matrix(coords, mode=mode, min_seq_separation=min_seq_separation, thresh=d_cutoff)
    length=np.shape(distances)[0]
    len_substructures=np.shape(substructures)[0]
    if length>len_substructures: #We are applying substructures from a smaller protein to analyze a larger protein, so we only keep the part of the larger protein that is encompassed by these substructures
    	distances=distances[0:len_substructures, 0:len_substructures]
    nsubstructures=np.shape(substructures)[2]
    scores=np.zeros((nsubstructures))
    for s in range(nsubstructures): 
        sub=substructures[:,:,s]
        participation=np.multiply(distances, sub)#gives the overlap between native substrucutres and this snapshot's contacts
        if mode=='distances':
        	scores[s]=np.mean(participation[np.nonzero(participation)])
        elif mode=='binary':
        	scores[s] = np.sum(participation)
        	
    return scores, distances

def Identify_nonnative_substructures(distances, native_contacts, filter_distance=2, d_cutoff=6, min_clustersize=7, contact_sep_thresh=7, min_to_comp_substructures=12):
    #First we pad the native contacts matrix so that contacts right next to native contacts are not counted as non-natives
    Filter=cp.deepcopy(native_contacts)
    for d in range(-filter_distance, filter_distance+1):
        im1_to_add=np.roll(native_contacts, d, axis=1)
        if d<0:
            im1_to_add[:, d:]=0
        else:
            im1_to_add[:, 0:d]=0
        
        im2_to_add=np.roll(native_contacts, d, axis=0)
        if d<0:
            im2_to_add[d:,:]=0
        else:
            im2_to_add[0:d, :]=0
        Filter=Filter+im1_to_add + im2_to_add
        Filter[np.where(Filter)]=1
    nonnative_distances=np.multiply(distances, 1-Filter)
    nonnative_distances[np.where(nonnative_distances==0)]=np.inf
    nonnative_contacts=np.zeros(np.shape(nonnative_distances))
    nonnative_contacts[np.where( nonnative_distances<=d_cutoff)]=1
    
    N_nonnative_substructures=0
    if np.sum(nonnative_contacts)>=min_to_comp_substructures:
        unused, nonnative_substructures=ClusterSubstructures.Identify_native_substructures(None, d_cutoff=d_cutoff,min_clustersize=min_clustersize,  contact_sep_thresh=contact_sep_thresh, native_contacts=nonnative_contacts, plot=False, verbose=False)
        N_nonnative_substructures=np.shape(nonnative_substructures)[2]
    return N_nonnative_substructures
        


##### Start the analysis here#######

### First, load the substructure file if it is specified, otherwise, identify the substructures de novo! ###

if Substructure_file!=None:
	unused1, unused2, unused3, native_distances, native_substructures = joblib.load(Substructure_file)
	native_contacts=np.zeros(np.shape(native_distances))
	native_contacts[np.where((native_distances<d_cutoff) & (native_distances!=0))]=1
	#note, by native_distances and native_contacts, we actually mean the distances/contacts that are found in the protein used to produce the substructure file
	#So for instance, if the substructure file is the full version of the protein of interest, then we define the contacts/distnaces in the full protein as the "native" ones
	#now, load the protein of interest's coordinates	
	coords, resis=Analyze_structures.read_PDB(native_file, 'CA')
	length=np.shape(coords)[0]
	if length < np.shape(native_substructures)[0]:
		native_substructures = native_substructures[0:length, 0:length, :]  #only keep whatever portions of substructures are acutally found in the current protein
		native_distances=native_distances[0:length, 0:length]  #same for native distances and contacts
		native_contacts=native_contacts[0:length, 0:length]
else:  #gotta identify those substructures de novo
    native_contacts, native_substructures=ClusterSubstructures.Identify_native_substructures(native_file, d_cutoff=d_cutoff, plot=False, min_clustersize=min_clustersize,contact_sep_thresh=contact_sep_thresh )
    coords, resis=Analyze_structures.read_PDB(native_file, 'CA')
    native_distances=Analyze_structures.compute_contacts_matrix(coords, mode='distances', min_seq_separation=8)


for N, directory in enumerate(directories):
    temps=temperatures[N]
    PDB_files=[]
    for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))     
    scores=np.zeros((len(PDB_files), np.shape(native_substructures)[2]))  #scores[f, n] tells you, for file f, what is the average distance between alpha carbons that were assigned to native substructure n    
    N_nonnative_substructures=np.zeros(len(PDB_files))  #now many nonnative substrucutres does file f have?
    print('{} files total'.format(len(PDB_files)))
    for f, file in enumerate(PDB_files):
        if np.mod(f,500)==0:
            print('{} files completed'.format(f))
        scores[f,:], distances = Score_snapshot(file, native_substructures, mode, d_cutoff)
        if ID_nonnatives:
        	N_nonnative_substructures[f]=Identify_nonnative_substructures(distances, native_contacts, filter_distance=2, d_cutoff=trap_d_cutoff, min_clustersize=min_clustersize,  contact_sep_thresh=contact_sep_thresh)
        else:
        	N_nonnative_substructures[f] = np.nan
        
        	
    if mode=='binary': #convert the native distances matrix into a contacts matrix
    	native_distances = Analyze_structures.compute_contacts_matrix(coords, mode='binary', min_seq_separation=8)
    
    joblib.dump([scores, N_nonnative_substructures, PDB_files, native_distances, native_substructures], open('{}/{}'.format(directory, output_filename), "wb"), compress=3 )