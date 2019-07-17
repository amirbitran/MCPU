#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 13:38:14 2018




@author: amirbitran
"""

import glob
import numpy as np
import joblib


directory='FABG/MultiUmbrella2'
log_files=glob.glob( '{}/*.log'.format(directory))


accepted=np.zeros(len(log_files))  #overall acceptance ratio for replica exchange
rejected=np.zeros(len(log_files))

n_successes=np.zeros((len(log_files), len(log_files)))  #on how many occasions has a node that started in state i ended in state j at the end of the 75 exchanges?
#n_attempts=np.zeros((len(log_files, len(log_files))))  #on how many occasions have nodes i and j attempted an exchagne? We kind of have to approximate this

for n, file in enumerate(log_files):
    #prev_accepted=0
    #prev_rejected=0
    for line in open(file, 'r').readlines():
        if line[0:6]=='myrank':
            fields=[x for x in line.split(' ') if len(x)>0]
            rank=int(fields[2].split(',')[0])
        if line[0:4]=='RPLC':                        
            line=line.rstrip('\n')
            fields=line.split(' ')
            fields=[f for f in fields if len(f)>0]

            #if accepted!=prev_accepted or rejected!=prev_rejected:  #this node participated in (at least one) exchange during the 75 exchanges...we treat this as half an exchange for each of the two partners
            if int(fields[1])>=25000000:  #start keeping track a quarter of way through once things have equilibrated
                #n_attempts[rank, rank+1]+=0.5
                #n_attempts[rank, rank+nodes_per_temp]+=0.5
                partner=int(fields[9].split('(')[0]) #with who's data did current node end up with after all 75 exchanges?
                if partner!=rank:
                    n_successes[rank, partner]+=1
            #prev_accepted=accepted
            #prev_rejected=rejected
           
            if int(fields[1])==99500000:  #at the end of the simulation, get overall acceptance/rejection rate for node
                accepted[rank]=int(fields[19].split(',')[0]) 
                rejected[rank]=int(fields[22].split(',')[0])
          
joblib.dump([n_successes,accepted, rejected], '{}/exchange_stats.dat'.format(directory))
