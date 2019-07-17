#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 08:13:37 2018

This script is intended to compute the number of native contacts for each snapshot for simulations where
this value is not already included in log file

@author: amirbitran
"""


import numpy as np
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import itertools
import joblib
import Analyze_structures
import sklearn
from sklearn import metrics
import LoopCluster
#from matplotlib import ticker
import glob
import pickle
import copy as cp
import natsort
import fnmatch
import copy as cp
import ClusterSubstructures
import visualize_PCA
#import Folding_rates
#import scipy.io



protein='FABG'
fileroot='fabg'


directories=['{}/replica_tight_spacing2'.format(protein)]
native_file='{}/files/{}_0.100_Emin.pdb'.format(protein, fileroot)
#native_file='{}/files/{}.pdb'.format(protein, fileroot)

temperatures=[['*.***']]

d_cutoff=6




def Compute_natives(snapshot, native_contacts, atom='CA', min_seq_separation=8 ):
    """
    Assigns a set of scores for a snapshot
    the ith score tells you what is the average distnace between pairs of residues residues that participate in the ith substructure, in this snapshto
    If the score is close to the characteristic contact distnace, then the substructure should be mostly formed
    """
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    contacts=Analyze_structures.compute_contacts_matrix(coords, min_seq_separation=min_seq_separation)
    return np.sum(np.multiply(contacts, native_contacts))





for N, directory in enumerate(directories):
    temps=temperatures[N]
    PDB_files=[]
    for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))     
    native_contacts, native_substructures=ClusterSubstructures.Identify_native_substructures(native_file, d_cutoff=d_cutoff, plot=False)
    natives=np.zeros(len(PDB_files))  #scores[f, n] tells you, for file f, what is the average distance between alpha carbons that were assigned to native substructure n    
    print('{} files total'.format(len(PDB_files)))
    for f, file in enumerate(PDB_files):
        if np.mod(f,500)==0:
            print('{} files completed'.format(f))
        natives[f]=Compute_natives(file, native_contacts)
    PDB_files, natives=visualize_PCA.sort_projections(PDB_files, natives)
    natives=natives.flatten()
    joblib.dump([natives, PDB_files], open('{}/N_natives.dat'.format(directory), "wb"), compress=3 )
