#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 13:34:04 2018
Randomly chooses snapshots with a given topological configuration from a simulation saved in directory, and 
prepares them as input files for an unfoldign simulation. They are stored in directory_name
@author: amirbitran
"""
import numpy as np
import os
import ClusterSubstructures
import ClusterPCA
import visualize_PCA
import shutil

protein='DHFR'
directory='MultiUmbrella2' #directory from which snapshots are drawn
temp='_0.8**_' #starting files are drawn from this temperature range
N_files=100
distinguish_traps = False
starting_config='1111100'
directory_name = 'starting_files_config_abcde'  #directory where select configurations are saved as 0.pdb, 1.pdb, etc..
thresh=1.7  #distance thresh for score assignment
min_step = 100000000 #only use steps that are greater than this number (inclusive)
max_step=np.inf  #Only use steps that are less than this number


def Get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times


scores,PDB_files=ClusterPCA.load_scores(protein, distinguish_traps=distinguish_traps, min_nonnative_substructures=2, directory=directory, thresh=thresh)



files_at_temp, s=visualize_PCA.Get_trajectory(scores, PDB_files, temp)
times = Get_times(files_at_temp)



#got rid of starting_cluster
candidates= [f for n,f in enumerate(files_at_temp) if s[n]==starting_config and times[n]>=min_step and times[n]<=max_step]

winners=np.random.choice(candidates, size=(N_files,))
print(winners)


starting_directory='{}/{}'.format(protein, directory_name)
if not os.path.isdir(starting_directory):
	os.mkdir(starting_directory)

for f, file in enumerate(winners):
	shutil.copy(file, '{}/{}.pdb'.format(starting_directory, f))
	
	
#for f, file in enumerate(winners):
#	shutil.copy('{}/{}.pdb'.format(starting_directory, f), file)
