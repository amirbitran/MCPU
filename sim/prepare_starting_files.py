# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:59:31 2017

@author: amirbitran

Randomly selects PDB files that contain a certain range of contacts,
and copies these into a directory that will contain starting structures for simulations
"""
import joblib
import numpy as np
import shutil 
import os 


print('Loading contacts data')
N=16   #Number of files to extract


directory='ADK'   #parent directory
PDB_prefix='adk'
contacts_min=900
contacts_max=1000
temp=0.65


initial_dir='{}/replica'.format(directory)

target_dir='{}/files/contacts_{}_{}'.format(directory, contacts_min, contacts_max)


if not os.path.exists(target_dir):
    os.mkdir(target_dir)

energies, contacts, temperatures=joblib.load("{}/All_data.dat".format(initial_dir))
temp_index= next(t for t in range(len(temperatures)) if temperatures[t]==temp)

pdb_indices=np.linspace(0,199000,200)


useful_indices=[]

for i in pdb_indices:
    if contacts_min<=contacts[int(temp_index),int(i)]<=contacts_max: useful_indices.append(int(i))


print('Selecting files')
chosen_indices=np.random.choice(useful_indices,N)   #randomly choose some PDB files


temp_str=str(temp)  #convert temp value into a string so that PDb file can be selected
while len(temp_str)<5:
    temp_str='{}0'.format(temp_str)  #add zeros to the end of temp string until there are 3 zeros following decimal


print('Copying selected files to directory')
for n,j in enumerate(chosen_indices):
    filename='{}_{}.{}'.format(PDB_prefix, temp_str, j*1000)
    #targetname='Contacts_{}.pdb'.format(contacts[int(temp_index),int(j)])
    targetname='{}.pdb'.format(int(n))
    #now just need to move the file
    shutil.copyfile('{}/{}'.format(initial_dir, filename),'{}/{}'.format(target_dir, targetname) )
        
