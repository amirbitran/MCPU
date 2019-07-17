#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 13:01:14 2018

FIND IDENTICAL TRAJECTORIES THAT ARE A PROBLEM!!!

@author: amirbitran
"""

import numpy as np
import visualize_PCA



def Find_identical_trajectories(data, PDB_files):
    #files, x=visualize_PCA.Get_trajectory(data, PDB_files, '*.***_**.')
    files, x=visualize_PCA.Get_trajectory(data, PDB_files, '*')
    x=x[:,0]
    
    
    x=np.reshape(x, (int(len(x)/200), 200))
    
    trajectories=[]
    for n in range(len(files)):
        if np.mod(n,200)==0:
            trajectories.append(files[n])
            
    
    id_trajectories=[]    
    
    for i in range(len(trajectories)):
        for j in range(i+1, len(trajectories)):
            if np.sum(x[i,:]-x[j,:])==0:
                id_trajectories.append((i,j))
    
    return id_trajectories


def Eliminate_identical_trajectories(data, PDB_files):
    """
    ASSUMES ALL TRAJECTORIES EVERYWHERE HAVE THE SAME LENGHT OF 200
    """
    id_trajectories=Find_identical_trajectories(data, PDB_files)
    for pair in id_trajectories:
        ind=pair[1]
        data[200*ind:200*ind+200,:]=np.nan*np.ones((200, np.shape(data) [1]))
        for i in range(200*ind, 200*ind+200):
            PDB_files[i]='NAN'
    
    data=data[~np.isnan(data).any(axis=1)]
    PDB_files=[f for f in PDB_files if f!='NAN']
    return data, PDB_files
        
        
        