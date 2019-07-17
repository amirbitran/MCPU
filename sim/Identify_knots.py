#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 08:13:37 2018

Used to find substructures with knots in a very heuristic way
 
Given the alpha carbon coordinates that comprise the loop, build a cube that contains the entire xyz domain encompassed by that loop
Then loop through all alpha carbons that comprise the thread. If any alpha carbons lives within this cube, then declare knot is formed

Of course this is a necessary condition but by no means sufficient...but the hope is that this, combined with substrucutre analysis, can be used to
declare whether a knot is present
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


loop_res = np.arange(88, 108, 1)  #which residues comprise the loop?
thread_min_res = 115  #all residues beyodn this one comprise the thread

directory = 'RSME/MultiUmbrella3'
output_filename = 'Knot_scores.dat'
fileroot = 'rsme'


temps = ['*.***']



PDB_files=[]

def Get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times
    
def Identify_knots(file,loop_res, thread_min_res):
    coords, resis=Analyze_structures.read_PDB(file, 'CA')
    loop = coords[loop_res,:]
    minx = np.min(loop[:,0])
    miny = np.min(loop[:,1])
    minz = np.min(loop[:,2])
    maxx = np.max(loop[:,0])
    maxy = np.max(loop[:,1])
    maxz = np.max(loop[:,2])        
    thread = coords[thread_min_res:, :] 
    knot = 0
    for ca in thread:
        if minx<ca[0]<maxx and miny<ca[1]<maxy and minz<ca[2]<maxz:
            knot = 1   
    return knot
            
for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))
    
times = Get_times(PDB_files)
         
scores=np.zeros(len(PDB_files))
print('{} files total'.format(len(PDB_files)))
for f, file in enumerate(PDB_files):
    if np.mod(f,500)==0:
        print('{} files completed'.format(f))
    scores[f] = Identify_knots(file, loop_res, thread_min_res)
        
    
joblib.dump([scores, PDB_files], open('{}/{}'.format(directory, output_filename), "wb"), compress=3 )



