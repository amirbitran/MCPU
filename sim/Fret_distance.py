#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 08:13:37 2018

Keeps track of distance between donor and acceptor amino acids
Need to fix this code so that you can just specify the amino acid numbers rather than the indexes
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



protein='HEMK_trunc11' #ex. MARR or MARR_trunc44
fileroot='hemk_t11' #ex. marr or marr_t44


directories=['{}/MultiUmbrella2'.format(protein)]  #what directories we analyze
temperatures=[['*.***']]
output_filename='Fret_distances.dat'



donor_atom='CE'
acceptor_atom='NZ'

donor_index = 0  #when we read the PDB for the given atom, we'll get a bunch of residues...which index corresponds to the donor? 0 for HEMK since donor is methionine
acceptor_index = 1  #Lysine 34 corresponds to the second (zero index 1) instance of CZ
#Note: Analyze_structures.read_PDB is written in a silly way so that a bunch of glycines get read when I try to read certain side chain atoms



for N, directory in enumerate(directories):
    temps=temperatures[N]
    PDB_files=[]
    for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))     
    print('{} files total'.format(len(PDB_files)))
    distances=np.zeros(len(PDB_files))
    for f, file in enumerate(PDB_files):
        if np.mod(f,500)==0:
            print('{} files completed'.format(f))
        donor_coords, donor_resis = Analyze_structures.read_PDB(file, donor_atom)
        acceptor_coords, acceptor_resis = Analyze_structures.read_PDB(file, acceptor_atom)
        donor_coord = donor_coords[donor_index]
        acceptor_coord = acceptor_coords[acceptor_index]
        distance = np.sqrt(np.sum((donor_coord - acceptor_coord)**2))     
        distances[f]=distance       
    joblib.dump([distances, PDB_files, donor_atom, acceptor_atom, donor_index, acceptor_index], open('{}/{}'.format(directory, output_filename), "wb"), compress=3 )
