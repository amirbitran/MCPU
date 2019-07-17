# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 17:59:20 2017

@author: amirbitran

Gets a pairwise distnace matrix between snapshots 
Various ways of defining pairwise distance

If mode='Contacts', then the distance between two snapshots is the total Manhattan distance between their contact matrices

if mode='TMscoring', then TM scores and RMSD values are computed

if mode='Dual', both modes are used.
"""

mode='Contacts'
#mode='tmscoring'
#mode='dual'
import Analyze_structures
import tmscoring  
import glob
import numpy as np
import pickle
#import timeit 


directory='ADK/unfolding'
PDB_files=glob.glob('{}/adk_0.850_*0000000'.format(directory))+glob.glob('{}/adk_0.850_*.0'.format(directory))   #all PDB files that are multiples of 10000
PDB_files.sort()
log=open('log.txt', 'w')

log.close()
#These tm scores are somewhat close than the ones given by the online server....although generally lower
#rmsd values differ from what pymol gives, although pymol may be aligning to optimize rmsd rather than TM score


def compute_pairwise_matrix(PDB_files, mode):
    """
    Gives you a pairwise matrix of TM scores--only computes lower triangle, since these scores are symmetric
    you input PDB_files, a list of the PDB files whose pairwise TM score you want to compute
    """
    TM=np.zeros((len(PDB_files), len(PDB_files)))
    RMSD=np.zeros((len(PDB_files), len(PDB_files)))
    
    delta_contacts=np.zeros((len(PDB_files), len(PDB_files)))
    
    #start = timeit.default_timer()
    
    
    
    for i in range(len(PDB_files)):
        TM[i,i]=1
        log=open('log.txt', 'a')
        for j in range(i):
        #for j in range(len(PDB_files)):  #computing everything for now since TM appears to be nonsymmetric
            
            log.write('Computing TM score and RMSD between files {} and {} \n'.format(i,j))
            P1=PDB_files[i]
            P2=PDB_files[j]
            if mode=='tmscoring':
            	sc = tmscoring.TMscoring(P1, P2) 
            	_, TM[i,j], RMSD[i,j] = sc.optimise()
            elif mode=='Contacts':
            	c1, resis=Analyze_structures.read_PDB(P1)  #first set of coordinates
            	c2, resis=Analyze_structures.read_PDB(P2)   #second set
            	
            	M1=Analyze_structures.compute_contacts_matrix(c1)   #contact matrix for first snapshot
            	M2=Analyze_structures.compute_contacts_matrix(c2)   #contact matrix for second snapshot
            	
            	delta_contacts[i,j]=np.sum(np.abs(M2-M1))
            elif mode=='dual':
            	sc = tmscoring.TMscoring(P1, P2) 
            	_, TM[i,j], RMSD[i,j] = sc.optimise()
            	c1, resis=Analyze_structures.read_PDB(P1)  #first set of coordinates
            	c2, resis=Analyze_structures.read_PDB(P2)   #second set
            	
            	M1=Analyze_structures.compute_contacts_matrix(c1)   #contact matrix for first snapshot
            	M2=Analyze_structures.compute_contacts_matrix(c2) 
            else: 
            	raise('Variable "mode" must be either, "tmscoring", "Contacts", or "dual"')
            	
            	
        log.close()
    return TM, RMSD, delta_contacts


TM, RMSD= compute_pairwise_matrix(PDB_files)
pickle.dump([PDB_files, TM, RMSD],open("{}/PDB_distances_{}".format(directory, mode), "wb"))





