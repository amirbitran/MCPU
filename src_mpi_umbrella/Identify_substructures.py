import sklearn
from sklearn import metrics
import LoopCluster
import glob
import pickle
import copy as cp
import numpy as np
#import visualize_PCA
#import natsort
#import fnmatch
import Analyze_structures
import copy as cp

files_path='/net/shakfs1/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/ADK/files'
structure='adk_0.100_Emin.pdb'

min_participants=4, 
atom='CA' 
d_cutoff=6, 
min_seq_separation=8

contact_sep_thresh=6
min_clustersize=6



coords, resis=Analyze_structures.read_PDB('{}/{}'.format(files_path, structure), atom)
   
native_contacts=Analyze_structures.compute_contacts_matrix(coords, thresh=d_cutoff, min_seq_separation=min_seq_separation)

positions=np.where(native_contacts==1) #which residues participate in contacts
positions=np.transpose(positions)

#print(positions)

M=metrics.pairwise.pairwise_distances(positions, metric='manhattan')  #how far is each contact from each other contact?

#To find connected components, I  use my loop_cluseter function by feeding in the positions ofr the contacts instead of the "files",
#as well as above matrix M as d_contacts

#If I want to be more fancy, I can probably find the connected ocmponents with a Laplacian matrix, but too lazy to do this now...
clusters, pairs_in_substructures, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(contact_sep_thresh, positions, M, sort_orphans=False, min_clustersize=min_clustersize)


#pairs in substructures is a list of sublists, each of which correspodns to a given substructure
#Within a given sublist, there are a bunch of pairs which tell you which pairs of residues belong to that substructure

#The above is in a messy form, so we convert it into a form that allows for numpy indexing,
#where we have a list of sublists, each sublist contains two arrays, the first of whcih gives the first indices for the interactin gresidues
pairs_in_substructures=[[np.array(C)[:,0], np.array(C)[:,1]] for C in pairs_in_substructures]



nsubstructures=len(pairs_in_substructures)  #we now produce a set of matrices...the ith page tells you which contacts belong to the ith substructure
substructures=np.zeros((np.shape(native_contacts)[0], np.shape(native_contacts)[1], nsubstructures))

substructures_file=open('{}/Substructures'.format(files_path), 'w')


for n in range(nsubstructures):
    SS=np.zeros(np.shape(native_contacts))
    nresidues=int(np.shape(SS)[0])
    SS[pairs_in_substructures[n]]=1
    SS= np.maximum( SS, SS.transpose() )  #make symmetric
    SS=SS.flatten()
    substructure_string=''
    for s in SS:
    	substructure_string='{}{}'.format(substructure_string, int(s))
    substructures_file.write('{}\n'.format(substructure_string))




substructures_file.close()	
