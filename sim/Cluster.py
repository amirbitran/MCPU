# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:27:32 2017

@author: amirbitran
"""

import numpy as np
import joblib
import matplotlib.pyplot as plt
import scipy
from scipy import linalg
import Analyze_structures
import pickle


#####################################################################################



def loop_cluster_contacts(thresh, files, d_contacts, sort_orphans=True):
    if d_contacts[1,0]!=d_contacts[0,1]:
        d_contacts=d_contacts+np.transpose(d_contacts)
    contacts_red=np.zeros(np.shape(d_contacts))
    for i, row in enumerate(d_contacts):
        for j, entry in enumerate(row):
            if entry<=thresh: contacts_red[i,j]=1
    
    
    #We start with the 0th structure. the goal is to find all nodes in the network that can be reached starting the 0th node entirely by traversing edges whose weight is greater than thresh
    
    #First, we find all nodes that only require one such passage
    
    
    clusters=[]
    clustered_indices=[]
    unsorted_indices=np.array(range(len(files)))
    
    
    
    while len(unsorted_indices)>0:
        curr_cluster=np.unique(np.nonzero(contacts_red[unsorted_indices[0],:])[0])
        
        
        
        new_layer=curr_cluster
        
        curr_len=len(curr_cluster)
        
        new_len=0
        
        
        while new_len != curr_len:
            curr_len=len(curr_cluster)
            for node in new_layer:
                new_layer=np.append(new_layer, np.nonzero(contacts_red[node,:])[0])
            new_layer=np.unique(new_layer)
            
            new_layer=np.array([x for x in new_layer if not np.any(curr_cluster==x)])
            
            curr_cluster=np.append(curr_cluster, new_layer)
            new_len=len(curr_cluster)
        
        
        curr_cluster=np.sort(curr_cluster)
        
        clusters.append(list(curr_cluster))
        clustered_indices.extend(curr_cluster)
        
        
        #clusters_zipped=[]
        unsorted_indices=np.array([x for x in range(len(files)) if x not in clustered_indices ])
    
    
    singlets=[x for x in clusters if len(x)==1]   #clusters with only one file--these are set aside
    clusters=[x for x in clusters if len(x)>1] #only keep the "fat" clusters with more than one file
    
    clusters_tally=np.zeros((len(files), len(clusters)))   #note the assignments in a matrix: rows are files, columns are fat clusters

    for n,x in enumerate(clusters):
        for y in x:
            clusters_tally[int(y),n]=1
    
    clustered_indices=np.where(np.sum(clusters_tally, axis=1)==1)[0]
    #unclustered_indices=np.where(np.sum(clusters_tally, axis=1)==0)[0]
    
    
    if sort_orphans:
        spots_for_singlets=np.zeros((len(files), len(clusters)))
        
        for n, singlet in enumerate(singlets):
            row=d_contacts[int(singlet[0]),:]
            indices=np.where(row==np.min(row[clustered_indices]))[0]   #yields all indices in row such that row(index)=minimum distance to current file among those files that have been classified 
            best_index=[i for i in indices if i in clustered_indices][0]    # we want
            cluster_assignment=np.where(clusters_tally[best_index,:]==1)[0][0]
        
            spots_for_singlets[int(singlet[0]), cluster_assignment]=1
        
        clusters_tally=clusters_tally+spots_for_singlets
        
    
    clusters=[]
    clustered_files=[]
    for n, col in enumerate(clusters_tally.T):
        clusters.append([i for i in range(len(col)) if col[i]!=0])
        clustered_files.append([files[i] for i in clusters[n]])
        
        print("Cluster {}: {} \n {} files total \n".format(n, clustered_files[n], sum(col)))
        
    #go through all pairs of files and figure out if two members of pair are in the same cluster or not.
    #If they are, add the distanec between these two files to the variable intracluster distances
    #Otherwise, add this distnce to intercluster_distances
    intracluster_distances=[]
    intercluster_distances=[]
    for i in range(len(d_contacts)):
        for j in range(i+1):
            if  np.size(np.where(clusters_tally[i,:])[0])!=0 and np.size(np.where(clusters_tally[j,:])[0])!=0:
                if np.where(clusters_tally[i])[0]==np.where(clusters_tally[j])[0]:  #they are in the same cluster
                    intracluster_distances.append(d_contacts[i,j])
                else:
                    intercluster_distances.append(d_contacts[i,j])
    
    mean_intracluster=np.mean(intracluster_distances)
    mean_intercluster=np.mean(intercluster_distances)
        
    print('Mean distance within clusters: {}'.format(mean_intracluster))
    print('Mean distance between clusters: {}'.format(mean_intercluster))
    

    return clusters, clustered_files, mean_intercluster, mean_intracluster
            
            
            
            
#####################################################################################


directory='ADK/unfolding'
temp='0.900'
files, d_contacts=joblib.load('{}/PDB_distances_{}_Contacts.dat'.format(directory, temp))
d_contacts=d_contacts+np.transpose(d_contacts)
thresh=68  #seems like may be a good value
clusters, clustered_files, mean_intercluster, mean_intracluster = loop_cluster_contacts(thresh, files, d_contacts)
#thresh=155
#clusters, clustered_files, cluster_distances, mean_intracluster, mean_intercluster = loop_cluster_contacts(thresh, files, d_contacts, assign_orphans=False)
        
        
        
reference_coords, resis=Analyze_structures.read_PDB(clustered_files[0][0])
reference_matrix=Analyze_structures.compute_contacts_matrix(reference_coords)


#we will now make a matrix where the ith page corresponds to the mean contacts map of the ith cluster
cluster_contact_maps=np.zeros((len(reference_matrix), len(reference_matrix), len(clustered_files)))


for n,cluster in enumerate(clustered_files):
	print('Reading files in cluster{}'.format(n))
	local_maps=np.zeros((len(reference_matrix), len(reference_matrix), len(cluster)))   #set of contact maps for this particular cluster
	for i, file in enumerate(cluster):
		coords, resis=Analyze_structures.read_PDB(file)
		local_maps[:,:,i]=Analyze_structures.compute_contacts_matrix(coords)	
	cluster_contact_maps[:,:,n]=np.mean(local_maps, axis=2)
	

pickle.dump([clusters, clustered_files, mean_intercluster, mean_intracluster, cluster_contact_maps],open("{}/Clusters_{}.dat".format(directory, temp), "wb"))