#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 09:09:54 2018

@author: amirbitran
"""



import numpy as np
#import pymbar # for MBAR analysis
#from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
import glob


##############################################################################
# First, read the file
##############################################################################

directory='1igd/MultiUmbrella4'

log_files=glob.glob('{}/*.log'.format(directory))



def get_temp(filename):
	splitline=filename.split(sep='/')
	split2=splitline[2].split('_')
	#print(split2)
	while split2[1][0] not in '0123456789':
		del split2[1]
	temp=float(split2[1][0:5])
	print(temp)
	return temp



exchanges=[]
times=[]
natives=[]
ranks=[]
temperatures=[]
setpoints=[]

for f, file in enumerate(log_files):
    curr_exchanges=[]
    curr_natives=[]
    temperatures.append(get_temp(file))
    openfile=open(file)
    

    for line in openfile.readlines():
        line=line.rstrip('\n')
        if len(line)>0 and line[0:6]=='myrank':
            entries=line.split()
            ranks.append(int(entries[2].rstrip(',')))

        

            
            
        if len(line)>0 and line[0:4]=='RPLC':
            entries=line.split()
            From=int(entries[9].split('(')[0])
            curr_exchanges.append(From)
            curr_natives.append(int(entries[7]))
            if f==0: times.append(int(entries[1]))
            
        if len(line)>0 and line[0:5]=='STEP ':
            entries=line.split()
            if entries[1]=='0':
                setpoint=int(entries[-1])
                setpoints.append(int(entries[-1]))
            
            
    natives.append(curr_natives)
    exchanges.append(curr_exchanges)
            

natives=np.array(natives)
exchanges=np.array(exchanges)

joblib.dump([exchanges, times, natives, ranks, temperatures, setpoints],'{}/Exchanges.dat'.format(directory))
            
 
        
        