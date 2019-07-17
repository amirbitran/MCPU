import numpy as np
#import pymbar # for MBAR analysis
#from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
#import matplotlib.pyplot as plt
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
import visualize_PCA
import ClusterSubstructures
 
 
"""
Loops through all PDB files that you care about* and saves either distances between every pair of residues in each file,
or contacts (a binary variable)
provided that the distances are at most d_cutoff

*By this, we mean that the PDB files shoudl be at some temperature of interest and potentially topological configuration 
"""  


native_file='DHFR/files/dhfr_0.100_Emin.pdb'
directory='DHFR/MultiUmbrella2'
fileroot='dhfr'
temps=['0.850', '0.875', '0.900']
scores_of_interest=['0000000'] #only read PDB_files whose score is contained in this array
#scores_of_interest='All'


#mode='binary'
min_seq_separation=8
atom='CA'

d_cutoff=8 #cutoff for contact distance #6.5 for FABG, 7.8 for MarR, 8 for HemK and protein G (with exception of ac contact maps, which is 10)
d_thresh=1.7  #how close must average distance between residues in substructures be (relative to native distance)for that substructure to be considered formed? 
#The following parameters are unnecessary, since Read_contact_maps does not actually identify substructures:
#min_clustersize=7 #7 for FABG
#contact_sep_thresh=4 #7 for FABG
distinguish_traps = False

min_step = 150000000
max_step = 200000000
#time_of_interest='3********'  #set to '*0' if we care about all times

distance_map_name = 'Distance_maps_null.dat'  #how we call the saved distance maps
score_file = 'Substructure_scores.dat'



def Get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times
    
def load_scores(directory,score_file='Substructure_scores.dat', distinguish_traps=True, thresh=2, elim_ID=True, convert_to_binary=True, min_nonnative_substructures=2, gmix=False):
    """
    Thresh tells you how far residues in a substructure have to be relative to native distance, on average, for substructure to be considered formed
    if elim_ID is true, you eliminate identical trajectories
    """
    try: 
        Scores, N_nonnative_substructures, PDB_files, native_contacts, substructures= joblib.load('{}/{}'.format(directory, score_file))  
        a, Scores=visualize_PCA.sort_projections(PDB_files, Scores)
        a, N_nonnative_substructures=visualize_PCA.sort_projections(PDB_files, N_nonnative_substructures)
        PDB_files=a
    except ValueError:
        Scores, PDB_files, native_contacts, substructures=joblib.load('{}/{}'.format(directory, score_file))
        PDB_files, Scores=visualize_PCA.sort_projections(PDB_files, Scores) 
    mean_substructure_distances=[]
    for i in range(np.shape(substructures)[2]):
        x=np.multiply(substructures[:,:,i], native_contacts)
        mean_substructure_distances.append(np.mean(x[np.where(x)]))
    mean_substructure_distances=np.array(mean_substructure_distances)
    if convert_to_binary:
        Scores=ClusterSubstructures.Substructure_scores_to_barcodes(Scores/mean_substructure_distances, thresh=thresh, gmix=gmix)
        #Scores=ClusterSubstructures.Substructure_scores_to_barcodes(Scores, thresh=thresh, gmix=gmix)
    if distinguish_traps:
        Scores=Identify_nonnative_states(Scores, N_nonnative_substructures, min_nonnative_substructures=min_nonnative_substructures)
    return Scores, PDB_files

def Identify_nonnative_states(scores, N_nonnative_substructures, min_nonnative_substructures=2):
    new_scores=[]
    N0=len(scores[0])
    for n,s in enumerate(scores):
        if N_nonnative_substructures[n]>=min_nonnative_substructures:
            new_scores.append('{}*'.format(s))
        else:
            new_scores.append(s)
    return new_scores
        

native_coords, resis=Analyze_structures.read_PDB(native_file, atom)
#we get a contact map with a min seq separation larger than usual to avoid helices
   
native_distances=Analyze_structures.compute_contacts_matrix(native_coords, mode='distances',thresh=d_cutoff,min_seq_separation=min_seq_separation)
#TODO: Optimize contact map calculation for computational efficiency     


scores, score_files=load_scores(directory, thresh=d_thresh, min_nonnative_substructures=1, score_file = score_file, distinguish_traps = distinguish_traps )


 
 
PDB_files=[]
#for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*pdb'.format(directory, fileroot, temp)) #USE ON  Mac
for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))

times = Get_times(PDB_files)
PDB_files = [file for f, file in enumerate(PDB_files) if times[f]>=min_step and times[f]<max_step ]


files_of_interest=[]

filescores=[]  #score for each file that was read

if scores_of_interest!='All':  
	for k, f in enumerate(PDB_files):
		ind=score_files.index(f)
		if scores[ind] in scores_of_interest:
			files_of_interest.append(f)
			filescores.append(scores[ind])
else:
	files_of_interest=PDB_files
        


print('{} files total'.format(len(files_of_interest)))

distance_maps=np.zeros((np.shape(native_distances)[0], np.shape(native_distances)[1], len(files_of_interest)))
  
  

      
for n,snapshot in enumerate(files_of_interest):
    if np.mod(n,1000)==0:
        print(n)
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    distance_maps[:,:,n]=Analyze_structures.compute_contacts_matrix(coords, mode='distances',thresh=d_cutoff, min_seq_separation=min_seq_separation)
    index=PDB_files.index(snapshot)
    #TODO! SUBTRACT OFF NATIVE CONTACTS!!!


contacts = np.zeros(np.shape(distance_maps)) 
contacts[np.where((distance_maps<d_cutoff) & (distance_maps!=0))]=1  #page i gives the contacts matrix for the ith snapshot

#we staple a page at the end of this contacts array that gives the average distnace between residue i and j, looking only that those snapshots where i and j form a bona fide contact
page = np.multiply(distance_maps, contacts)
page[np.where(page==0)]=np.nan
page = np.nanmean(page, axis=2)


contacts= np.dstack((contacts, page ))
            
joblib.dump([contacts, files_of_interest, filescores ], open('{}/{}'.format(directory, distance_map_name), 'wb'), compress=6)
    



