#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:13:32 2018

This script loops through a bunch of snapshots of dimer simulations and counts the
number of contacts at the interface, as well as the number of substructures
involving interface residues

@author: amirbitran
"""

import numpy as np
import joblib
import Analyze_structures
import glob
import ClusterSubstructures

mon1=np.arange(0, 144, 1) #residue indices (zero indexing) for monomer 1
mon2 = np.arange(158,302, 1) #monomer 2

directory='MARR_dimer/MultiUmbrella2'
fileroot='marr_d'
temperatures = ['*.***']

output_filename = 'Interface_contacts.dat'

d_cutoff = 8.5  #cutoff for substrucutres...since these are helical contacts, we are lenient
min_clustersize=7 #7 for MARR, 7 for FABG
contact_sep_thresh=3 #7 for FABG, 3 for MARR
min_seq_separation=8

PDB_files = []
for temp in temperatures: PDB_files += glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))  


n_contacts = np.zeros(len(PDB_files)) #number of interface contacts
n_subs = np.zeros(len(PDB_files)) #number of interface substructures
print('{} files total'.format(len(PDB_files)))

#loop through snapshots
for n, snapshot in enumerate(PDB_files):
    if np.mod(n,500)==0:
        print('{} files completed'.format(n))
    coords, resis=Analyze_structures.read_PDB(snapshot, 'CA')
    distances=Analyze_structures.compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
    
    x = distances[mon2, :]
    x = x[:, mon1]  #at this point, x is the distance map for pairs of residues (i,j) where i is in monomer 2 and j is in monomer 1
    
    contacts = np.zeros(np.shape(x))
    contacts[np.where(x <= d_cutoff)]=1
    n_contacts[n] = np.sum(np.sum(contacts))
    if n_contacts[n]>=min_clustersize: #there's a CHANCE of a substructure. otherwise, no chance!
    	unused, subs=ClusterSubstructures.Identify_native_substructures(None, d_cutoff=d_cutoff,min_clustersize=min_clustersize,  contact_sep_thresh=contact_sep_thresh, native_contacts=contacts, plot=False, verbose=False)
    	n_subs[n] = np.shape(subs)[2]
    else:
    	n_subs[n] = 0
    	

joblib.dump([n_contacts, n_subs, PDB_files], open('{}/{}'.format(directory, output_filename), "wb"), compress=3 )