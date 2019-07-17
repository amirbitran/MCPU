# -*- coding: utf-8 -*-
"""
Created on Wed May  2 17:53:18 2018

@author: amirbitran
"""

import numpy as np
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
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
import natsort
import fnmatch
import ClusterPCA
import copy as cp
import scipy.io


def Identify_native_substructures(native_file, d_cutoff=6, min_seq_separation=8, atom='CA', contact_sep_thresh=7, min_clustersize=6, plot=True, native_contacts=[], verbose=True):
    """
    Identify substructures within native file contact map
    Using the following strategy
    
    We produce a contact map which is a bunch of dots
    Contacts correspond to pairs of residues that are less than d_cutoff apart
    6 Angstroms is generally a good value, but may want a larger value for helical proteins where residues interact primarily
    via sidechains, and thus the alpha carbons are further apart
    
    We only count contacts if the residues are separated by min_seq_separation along the primary sequence
    We set min_seq_separation relatively high because we don't care to include intra-helix contacts within our contact map
    
    Ce can calculate the "Manhattan" distnace between every pair of dots on that contact map 
    and build a graph of contacts in which two contacts are connected if they are less than some
    threshold distance, contact_sep_thresh, apart in the contact map 
    
    Then, we find all connected components of this graph, each of which is a substructure
    But we only keep substructures whose number of residues is at least min_clustersize, to avoid sparse contacts here and there that we dont' care about
    
    Gives you option to input native contacts a priori, but by defualt you don't do this (value is set to None)
    """
    if len(native_contacts)==0:
        coords, resis=Analyze_structures.read_PDB(native_file, atom)
        #we get a contact map with a min seq separation larger than usual to avoid helices
   
        native_contacts=Analyze_structures.compute_contacts_matrix(coords, thresh=d_cutoff, min_seq_separation=min_seq_separation)
        #TODO: Optimize contact map calculation for computational efficiency 
    
    #Plot the native contact map
    #if plot:
    #    plt.figure()
    #    plt.imshow(native_contacts, cmap='Greys') 
    #    plt.title('Native contacts for {}'.format(native_file))
    
    positions=np.where(native_contacts==1) #which residues participate in contacts
    positions=np.transpose(positions)

    M=metrics.pairwise.pairwise_distances(positions, metric='manhattan')  #how far is each contact from each other contact?
    
    #To find connected components, I  use my loop_cluseter function by feeding in the positions ofr the contacts instead of the "files",
    #as well as above matrix M as d_contacts
    
    #If I want to be more fancy, I can probably find the connected ocmponents with a Laplacian matrix, but too lazy to do this now...
    clusters, pairs_in_substructures, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(contact_sep_thresh, positions, M, sort_orphans=False, min_clustersize=min_clustersize, verbose=verbose)
    

    #pairs in substructures is a list of sublists, each of which correspodns to a given substructure
    #Within a given sublist, there are a bunch of pairs which tell you which pairs of residues belong to that substructure
    
    #The above is in a messy form, so we convert it into a form that allows for numpy indexing,
    #where we have a list of sublists, each sublist contains two arrays, the first of whcih gives the first indices for the interactin gresidues
    pairs_in_substructures=[[np.array(C)[:,0], np.array(C)[:,1]] for C in pairs_in_substructures]
    
    
    
    nsubstructures=len(pairs_in_substructures)  #we now produce a set of matrices...the ith page tells you which contacts belong to the ith substructure
    substructures=np.zeros((np.shape(native_contacts)[0], np.shape(native_contacts)[1], nsubstructures))
    for n in range(nsubstructures):
        SS=np.zeros(np.shape(native_contacts))
        SS[pairs_in_substructures[n]]=1
        substructures[:,:,n]=SS
    if plot:
        Visualize_substructures(native_contacts, substructures)
    #print(positions)
    return native_contacts, substructures





#Here's the deal: This snippet of code computes a "score" for each substructure--that is, the average distance between the apirs of residues participating in that substructure
#snapshot='FABG_tighter_temps/fabg_0.570.0.pdb'
#snapshot='FABG_tighter_temps/fabg_0.660.94500000.pdb'




def Assign_snapshot(snapshot, substructures, min_participants=4, atom='CA', d_cutoff=6, min_seq_separation=8):
    """
    Identifies which substructures are present in a given snapshot (ex. 'adk_0.810_0.pdb')
    
    To do this, we loop through all pages of the substructure array, each of 
    which contains 1 at contact map positions participating in a given substructure
    
    We then take the dot product of that page with the contact map for the snapshot
    to find the number of contacts in that snapshot belonging to that substructure
    
    If that number of contacts is at least min_participants, then we say the substructure is present in the snapshot
    
    Returns a string, ex "10010" which is 1 at indicies correpsonding substructures that are present in the sanpshot,
    and 0 at other positions
    
    """
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    
    #TODO: Optimize the following step for computational efficiency
    contacts=Analyze_structures.compute_contacts_matrix(coords, thresh=d_cutoff, min_seq_separation=min_seq_separation)
    nsubstructures=np.shape(substructures)[2]
    Assignment=''
    participation=np.zeros(np.shape(substructures))
    n_participants=np.zeros((nsubstructures))
    for s in range(nsubstructures): 
        sub=substructures[:,:,s]
        participation[:,:,s]=np.multiply(contacts, sub)#gives the overlap between native substrucutres and this snapshot's contacts
        n_participants[s]=np.sum(participation[:,:,s])
        if n_participants[s]>=min_participants:
            Assignment='{}1'.format(Assignment)
        else:
            Assignment='{}0'.format(Assignment)


    return Assignment, n_participants






def Read_contact_maps(native_file, directories, temperatures, fileroot, mode='binary',thresh_distance=6,min_seq_separation=8, atom='CA'):
    """
    Loops through all PDB files that you care about and saves either distances between every pair of residues in each file,
    or contacts (a binary variable)
    provided that the distances are at most thresh_distance
    
    Never ended up using this since it involves too much memory

    """
    native_coords, resis=Analyze_structures.read_PDB(native_file, atom)
    #we get a contact map with a min seq separation larger than usual to avoid helices
   
    native_distances=Analyze_structures.compute_contacts_matrix(native_coords, mode=mode,thresh=thresh_distance,min_seq_separation=min_seq_separation)
    #TODO: Optimize contact map calculation for computational efficiency     
    
    
    PDB_files=[] #first, we figure out what PDB files we need
    for N, directory in enumerate(directories):
        temps=temperatures[N]
        if  os.path.exists('{}/Distance_maps.dat'.format(directory)):
            [distance_maps, PDB_files, native_distances]=joblib.load('{}/Distance_maps.dat'.format(directory))   
            print('Loaded data from directory {}'.format(directory))
            		
            original_distance_maps=cp.deepcopy(distance_maps)  #copy this variable since we will need it later in its intact form
            original_files=PDB_files
            		
            #First thing: We see if there are any files and maps in here that we don't need
            #needed_files are all the files we need
            needed_files=[]
            #for temp in temps: needed_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp))  USE IN ODYSSEY
            for temp in temps: needed_files+=glob.glob('{}/{}_{}*.*pdb'.format(directory, fileroot, temp)) 
            				
            indices_to_remove=[i for i in range(len(PDB_files)) if PDB_files[i] not in needed_files]
            		
            		
            distance_maps=np.delete(distance_maps, indices_to_remove, 2) #removes the pages (hence along axis "2") indicated by the above indices
            PDB_files = [PDB_files[i] for i in range(len(PDB_files)) if i not in indices_to_remove]
            		
            		    	
            #We need to double check that all the files we want are in Distance_maps.dat
            #If some are missing, we should add them using Get_aligned_files function
    
            #Now check if all files in PDB files are in needed_files
            #missing_files is a list of files that are in needed_files but not in PDB_files
            		
            missing_files=[f for f in needed_files if f not in PDB_files]
            		
            #Get the missing distance maps from missing_files, and append accordingly
            		
            if len(missing_files)>0:
                print('Some files are missing!')
                missing_data=np.zeros((np.shape(native_distances)[0], np.shape(native_distances)[1], len(missing_files)))
                for n, snapshot in enumerate(missing_files):
                    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
                    missing_data[:,:,n]=Analyze_structures.compute_contacts_matrix(coords, mode=mode,thresh=thresh_distance, min_seq_separation=min_seq_separation)
                
                print('Obtained data for additional PDB files')

                Data_to_save=np.concatenate((original_distance_maps, missing_data), axis=2)  #we combine the original data (without having removed any files) with the data we accumulated during this run, and save everything, so that our repository grows
                files_to_save=original_files + missing_files
                joblib.dump([Data_to_save, files_to_save, native_distances], open("{}/Distance_maps.dat".format(directory),"wb"), compress=6)
            else:
                joblib.dump([distance_maps, PDB_files, native_distances], open('{}/Distance_maps.dat'.format(directory), 'wb'), compress=6)
            			
		
        else: 
            for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*pdb'.format(directory, fileroot, temp)) #USE ON  Mac
            #for temp in temperatures: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) #USE ON  Odyssey

            distance_maps=np.zeros((np.shape(native_distances)[0], np.shape(native_distances)[1], len(PDB_files)))
            
            for n,snapshot in enumerate(PDB_files):
                coords, resis=Analyze_structures.read_PDB(snapshot, atom)
                distance_maps[:,:,n]=Analyze_structures.compute_contacts_matrix(coords, mode=mode,thresh=thresh_distance, min_seq_separation=min_seq_separation)
                
            joblib.dump([distance_maps, PDB_files, native_distances], open('{}/Distance_maps.dat'.format(directory), 'wb'), compress=6)
        
    


def Visualize_substructures( native_contacts, substructures):
    """
    Visualizes substructures as follows
    Everything that is a native contact but not part of any substructure will have value -1 on shown image
    Meanwhile, everything that is part of substructure i (i ranges from 0 to N_substructures-1) will have value i
    Finally, all non-contacts will just be Nans and appear white
    """
    plt.figure()
    #Cmap_masked=np.ma.masked_array(native_contacts, pairs_in_substructures[s])
    substructure_image=np.zeros(np.shape(native_contacts))
    for s in range(np.shape(substructures)[2]):      
        substructure_image+=(s+1)*substructures[:,:,s]
        

        #substructure_image[pairs_in_substructures[s]]=s+1
        #im1=plt.imshow(template,alpha=0.7, cmap='jet', interpolation='none')
    im=substructure_image+native_contacts-2   #ensures unassigned contacts have value -1, while assign contacts have value corresponding to their substructure number
    im[im==-2]=np.nan
    plt.imshow(im, cmap='jet')
    plt.tick_params(labelsize=30)
    
    
    colors=cm.get_cmap('jet')
    
    for s in range(np.shape(substructures)[2]):      
        #Let's annotate
        y_pos=np.where(substructures[:,:,s])[0][0]-3
        x_pos=np.where(substructures[:,:,s])[1][0]+3
        curr_color=colors((s+1)/(np.max(substructure_image) ))
        plt.annotate('{}'.format(s), (x_pos, y_pos), fontsize=30, color=curr_color)

    
    nsubstructures=np.shape(substructures)[2]
    nbins=nsubstructures+1 #number of colors we are showing...add 1 to accoutn for unassigned contacts


    #The following if you need a colorbar
    #cb=plt.colorbar()
    #tick_locator = ticker.MaxNLocator(nbins=nbins)
    #cb.locator = tick_locator
    #cb.update_ticks()
    #plt.show()


def Visualize_substructure(s, native_contacts, pairs_in_substructures):
    """
    Visualize just a single substrucutre
    a bit messy, this hasn't been updated
    """
    plt.figure()
    #Cmap_masked=np.ma.masked_array(native_contacts, pairs_in_substructures[s])
    template=np.zeros(np.shape(native_contacts))
    template[pairs_in_substructures[s]]=1
    #im1=plt.imshow(template,alpha=0.7, cmap='jet', interpolation='none')
    im1=plt.imshow(template+native_contacts, cmap='jet', interpolation='none')
    #im2=plt.imshow(native_contacts, cmap='Greys', interpolation='none')
    plt.show()
    
    
    



def Assign():   
    native_file='ADK/replica_tight_spacing3/adk_0.600.0'
    directories=['ADK/replica_tight_spacing3', 'ADK/unfolding']
    temperatures=[['*.***'], ['*.***']]
    fileroot='adk'
    
    print('Obtaining native substructures')
    native_contacts, substructures=Identify_native_substructures(native_file, plot=False)
    PDB_files=[]
    print('Obtaining list of PDB files')
    for N, directory in enumerate(directories):
        temps=temperatures[N]
        for temp in temps: PDB_files+=glob.glob('{}/{}_{}*.*0'.format(directory, fileroot, temp)) #USE ON  Odyssey
    print('Obtaining assignments')
    Assignments=[]
    for n,snapshot in enumerate(PDB_files):
        Assignments.append(Assign_snapshot(snapshot, substructures))
        if np.mod(n,500)==0: print('{} files completed'.format(n))
    joblib.dump([Assignments, PDB_files, native_contacts, substructures ],open('{}/Substructure_assignments.dat'.format(directories[0]), 'wb'))
        
#     
#def lookup_files(label,labels, PDB_files, traj='All'):
#    """
#    returns all files that are assigned to a given label   
#    You can use traj to specify that you only care about files in a given trajectory
#    The defalt is 'All' (all trajectories are considered )
#    """
#    if np.isnan(label):
#        files=[f for n,f in enumerate(PDB_files) if np.isnan(labels[n])]
#    else:
#        files=[f for n,f in enumerate(PDB_files) if labels[n]==label ]
#    if traj!='All':
#         files = [f for f in files if fnmatch.fnmatch( f, '*_{}*'.format(traj))]
#    return files        

 ##############################
 
 
 
def Hub_clustering(protein,req_frac=0.96, combine_thresh=1000 ):
    """
    Kinetic hub-based clustering: A proposed algorithm, based on the observation that many unfolding trajectories pass through a relatively small set of configurations—these overrepresented configurations will be called “hubs”
    
    1. Filter oscillations in all trajectories to get rid of as many extraneous transitions as possible
    2. Construct some sort of kinetic distance matrix between every pair of points (simple transition matrix, mean waiting time matrix, etc)
    3. Rank the configurations by how often they occur in all unfolding trajectories
    
    Now begin the iterative algorithm
    
    1. Choose the most overrepresented configurations. Then find all configurations that are within a specified “distance” from this hub (different ways of defining this)…these configurations will be assigned to this hub
    
    2. Check if at a fraction of all points given by req_frac has been assigned to a cluster . If so, end the algorithm
    
    3. If not, add the next most represented hub, and repeat step 1, assigning each configuration to its NEAREST hub in case certain configurations are within the cutoff distance of multiple hubs
    
    4. Repeat step 3 until 2 is satisfied
    
    Then finally, we combine clusters whose hubs are within combine_thresh of each other
    """
    try:
        unfolding_points, PDB_files=ClusterPCA.load_scores(protein, directory='unfolding', convert_to_binary=False)
        times=np.array(ClusterPCA.Get_times(PDB_files))
        labels, unique_points=Bernie_elimination(unfolding_points, PDB_files, 0, plot=False)
    except TypeError:
        unfolding_points, PDB_files=ClusterPCA.load_scores(protein, directory='unfolding', convert_to_binary=True)
        times=np.array(ClusterPCA.Get_times(PDB_files))
        labels, unique_points=Bernie_elimination(unfolding_points, PDB_files, 0, plot=False)
        
    labels=np.array(labels)
    
    
    ###### filter trajectories #############
    
    zero_times=np.where(times==0)[0]
    
    for i, t in enumerate(zero_times):
        if t==zero_times[-1]:
            labels[t:]=ClusterPCA.filter_oscillations(labels[t:], thresh=7)
        else:
            nextt=zero_times[i+1]
            labels[t:nextt]=ClusterPCA.filter_oscillations(labels[t:nextt], thresh=7)
   
    #redo Bernie elimination now that we've potetnially gtten rid of some points
    unfolding_points=[unique_points[l] for l in labels]
    labels, unique_points=Bernie_elimination(unfolding_points, PDB_files, 0, plot=False) 
    unique_labels=np.unique(labels)
    labels=np.array(labels)
        
    ###### Sort hubs by representation #############
    
    rep=np.array([len(np.where(labels==u)[0]) for u in unique_labels])
    sorted_indices=np.argsort(-rep)
    sorted_hubs=unique_labels[sorted_indices]   
    
    
    ########### Construct a distance matrix #####################
        
    
    n=len(unique_points)
    #build transition matrix...am doing this myself rather than using the Markov_transmat function in Folding_rates, since that function is UGLY
    counts=np.zeros((n,n))
    norms=np.zeros((n,n))
    
    for t, i in enumerate(labels):
        if t<len(labels)-1 and times[t+1]!=0:
            i=int(i)
            j=int(labels[t+1])
            counts[i,j]+=1
            norms[i,:]+=1
    transmat=np.divide(counts, norms)
           
    
    #now we compute that distance metric
    
    labels=np.array(labels)
    distance=np.nan*np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            P_i_given_ij=len(np.where(labels==i)[0])/(len(np.where(labels==i)[0])+len(np.where(labels==j)[0]))
            P_j_given_ij=1-P_i_given_ij
            
            distance[i,j]=  P_i_given_ij/transmat[i,j] + P_j_given_ij/transmat[j,i]
            #distance[i,j]=  (transmat[i,j]+transmat[i,i]) *P_i_given_ij/transmat[i,j] + (transmat[j,i]+transmat[j,j]) *P_j_given_ij/transmat[j,i]
            
    distance[np.diag_indices(n)]=0  #sometimes it isn't small...
    
    
    ####### Begin the iterative algorithm--keep going until 95% of timepoints have been assigned, and the native state has been assigned #########
    
    #go=True
    N=0  #how many hubs to use in current iteration
    
    assigned=[] #labels that have been assigned to a cluster already already
    frac_assigned=0
    
    #while frac_assigned<0.95 or 0 not in assigned:
    while frac_assigned<req_frac:
        assigned=[]
        N+=1
        hubs=sorted_hubs[0:N]
        clusters=[[] for j in range(N)]
        
        for label in unique_labels:
            distances_to_hubs=distance[label, hubs]
            if np.min(distances_to_hubs)<np.inf:
                assignment=np.argmin(distances_to_hubs)
                clusters[assignment].append(label)
                assigned.append(label)

        frac_assigned=np.sum(rep[assigned])/np.sum(rep)
        
    
    #####This is something new I am attempting on Monday 8/13...namely doing loop clustering just on the kinetic hubs ####
    
    hub_distances=np.array([[d for d in distance[h,hubs]]for h in hubs])
    
    unused_output, who_to_merge, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(combine_thresh, np.arange(0,N,1), hub_distances, sort_orphans=False, min_clustersize=1)  
        
    new_clusters=[] #clusters after mergers
    
    #First, loop through clusters that have been merged
    for w in who_to_merge:
        new_cluster=[]
        for clust in w:
            for conf in clusters[clust]:
                new_cluster.append(conf)
        new_clusters.append( list(np.sort(new_cluster) ))
    
    #now loop through the clusters that have not been merged
    flat_list=[c for  w in who_to_merge for c in w]
    unmerged=[c for c in np.arange(0,N,1) if c not in flat_list]
    for u in unmerged:
        new_clusters.append(clusters[u])

    clusters=new_clusters
    
    #Sort cluster order by mean number of substructures
    n_substructures=[]
    for i,c in enumerate(clusters):
        n_substructures.append([])
        for l in c:
            n_substructures[i].append(np.sum([int(s) for s in unique_points[l]]))
            
    mean_n=np.array([np.mean(n) for n in n_substructures])
    order=np.argsort(-1*mean_n)
    clusters=[clusters[o] for o in order]
    N=len(clusters)
    
    
    #to which cluster does each substructure belong? If none, we assign a nan
    key=[] #the ith element tells you to which cluster do points with label i belong
    for p in unique_labels:
        if p in assigned:
            key.append([c for c in range(N) if p in clusters[c]][0])
        else:
            key.append(np.nan)
    
    unfolding_clusters=[ unique_points.index(p) for p in unfolding_points] #convert to numerical labels
    unfolding_clusters=np.array([key[p] for p in unfolding_clusters]) #convert to cluster assignments
    
    
    #Filter oscillations and get rid of Nans (unassigned points) by assigning them to the cluster that came immediately before
    
    for i, t in enumerate(zero_times):
        if t==zero_times[-1]:
            unfolding_clusters[t:]=ClusterPCA.filter_nans(unfolding_clusters[t:])
            unfolding_clusters[t:]=ClusterPCA.filter_oscillations(unfolding_clusters[t:], thresh=7)
        else:
            nextt=zero_times[i+1]
            unfolding_clusters[t:nextt]=ClusterPCA.filter_nans(unfolding_clusters[t:nextt])
            unfolding_clusters[t:nextt]=ClusterPCA.filter_oscillations(unfolding_clusters[t:nextt], thresh=7)
    
    
    ##########################End 8/13 strategy ###########
    
    clusters=[[unique_points[l] for l in clust] for clust in clusters]
    
    return unfolding_clusters, clusters, key, PDB_files
    
    




   #We've assigned most snapshots to clusters...now we need to asign the rest. To do this, we loop throuhg the unassigned snapshots, and 
   #assign each to whichever cluster that contains an assigned snapshot that is closest kinetically to the unassigned snapshot in question
#   
#
#    ######### A method I ended up not using #########
#
#
#    #My worry about this approach: Suppose there is some pathological scenario like the following:
#    #Suppose confiigs 0 and 3 have not been assigned...Suppose for whatever reason, the only cluster to which 0 has a connection is cluster 2 (a mostly unfolded cluster) maybe like a connction of 100000 or something (due to some accidental fluctuation)
#    #Suppose, though that 0 and 3 are closely connected, and that config 3 does have a connection to the folded cluster to which it shoudl belong (call it cluster 1)
#    #Then if 0 participates in part of the algorithm before 3, then 0 will get assigned (abbaranetly) to cluster 2, while 3 will get assigned (rightfully to cluster 1)
#    #But if we were to now switch so that 3 is assigned before 0, then both will rightfully end up in cluster 1
#    
#    #The solution may be to do this iteratively without actually adding the orphans to the clusters until the algorithm is complete
#   run=True
#   unassigned=[a for a in unique_labels if a not in assigned]
#   
#
#   while run:
#       for s in unassigned:
#           assigned_potential_partners=np.array([p for p in assigned if distance[s,p]<np.inf])  #is there anybody who has already been assigned to a cluster that 
#           chosen_partner=assigned_potential_partners[np.argmin(distance[s, assigned_potential_partners])]
#       
#    
   
   
#   
##    
#    ### this is the strategy for combining clusters that I was using on Friday 8/10 #####
#    proximity_thresh=300
#    overlap=np.zeros((N,N))
#    #For each pair of clusters, ask what fraction of configurations in a given cluster is within proximity_thresh of the partner cluster's hub
#    for i in range(N):
#        for j in range(i+1, N):
#            ntotal=len(clusters[i])+len(clusters[j])
#            n_proximal=0
#            for l in clusters[i]:
#                if distance[l, hubs[j]]<proximity_thresh:
#                    n_proximal+=1
#            for l in clusters[j]:
#                if distance[l, hubs[i]]<proximity_thresh:
#                    n_proximal+=1
#            overlap[i,j]=n_proximal/ntotal
#            #overlap[j,i]=overlap[i,j]
#            
#    
#    #Another method to combine clusters
##    proximity_thresh=100 #for each pair of clusters, count how many pairs of labels are within this kinetic distance of each other
##    overlap=np.zeros((N,N))
##    
##    for i in range(N):
##        for j in range(i+1, N):
##            n_pairs=len(clusters[i])*len(clusters[j])
##            number_proximal=0
##            for k in clusters[i]:
##                for l in clusters[j]:
##                    if distance[k,l]<proximity_thresh:
##                        number_proximal+=1
##            overlap[i,j]=number_proximal/n_pairs
##            overlap[j,i]=overlap[i,j]
#    
#    overlap_thresh=0.25 #to be combined, a pair of clusters must overlap in at least this fraction of pairs of points
#    overlap[np.diag_indices(N)]=1
#
#    unused_output, who_to_merge, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(1-overlap_thresh, np.arange(0,N,1), 1-overlap, sort_orphans=False, min_clustersize=1)   
#    
#    new_clusters=[] #clusters after mergers
#    
#    #First, loop through clusters that have been merged
#    for w in who_to_merge:
#        new_cluster=[]
#        for clust in w:
#            for conf in clusters[clust]:
#                new_cluster.append(conf)
#        new_clusters.append( list(np.sort(new_cluster) ))
#    
#    #now loop through the clusters that have not been merged
#    unmerged=[c for c in np.arange(0,N,1) if c not in np.array(who_to_merge).flatten() ]
#    for u in unmerged:
#        new_clusters.append(clusters[u])
#
#    clusters=new_clusters
#    
#    #Sort cluster order by mean number of substructures
#    n_substructures=[]
#    for i,c in enumerate(clusters):
#        n_substructures.append([])
#        for l in c:
#            n_substructures[i].append(np.sum([int(s) for s in unique_points[l]]))
#            
#    mean_n=np.array([np.mean(n) for n in n_substructures])
#    order=np.argsort(-1*mean_n)
#    clusters=[clusters[o] for o in order]
#    N=len(clusters)
#    
#    
#    #to which cluster does each substructure belong? If none, we assign a nan
#    key=[] #the ith element tells you to which cluster do points with label i belong
#    for p in unique_labels:
#        if p in assigned:
#            key.append([c for c in range(N) if p in clusters[c]][0])
#        else:
#            key.append(np.nan)
    
#    
#    #np.array(who_to_merge).flatten()
#    
#    
#    new_clusters=[]
#    have_I_been_combined=[False for c in clusters]
#    I_was_combined_with=[np.nan for c in clusters]
#    import copy as cp
#    for i,c in enumerate(clusters):
#        if not have_I_been_combined[i]:
#            have_I_been_combined[i]=True
#            combine_with=np.where(overlap[i,:]>=overlap_thresh)[0]
#            new_cluster=cp.deepcopy(c)
#            
#            
#            for z in combine_with:
#                if not have_I_been_combined[z]:
#                    for s in clusters[z]:
#                        new_cluster.append(s)
#                    have_I_been_combined[z]=True
#                else: #somebody else already incorporated cluster z
#            new_clusters.append(new_cluster)
#            have_I_been_combined[i]=True
#    
#    clusters=new_clusters
#    N=len(clusters)
#        #for z in combine_with:
#        #    clusters.remove(clusters[z])
            
        

        
    
def Kinetic_clustering(protein, thresh ):
    """
    Clusters topological configurations based on whether they are kinetically connected to each other
    Here's how we do it
    
    First we load all substructure scores and convert to barcode format
    Then we build a Markov tarnsition matrix for the transitions between the current scores
    Then, we construct a distance between configurations i and j by asking, given that you are in i or j, how long do you have to wait, on average, to transition to the other in the pair?
    Finally we use the threshold thresh to do a loop clustering, to find clusters of configurations that are all connected to each other by edges of distance less than thresh
    """

    unfolding_points, PDB_files=ClusterPCA.load_scores(protein, directory='unfolding', convert_to_binary=True)
    
    #if protein=='ADK_trunc26':  #only use higher temperatures, since lower temps have too many oscillations..
        #PDB_files, unfolding_points=visualize_PCA.Get_trajectory(unfolding_points, PDB_files, ['0.8**_', '0.9**_'])
    
    labels, unique_points=Bernie_elimination(unfolding_points, PDB_files, 0, plot=False)        
    
    
    
    x=[] #how many substrucutres are formed at each time
    for p in unfolding_points:
        x.append(np.sum([int(i) for i in p]))
    
    #make a .mat file...
    times=ClusterPCA.Get_times(PDB_files)
    scipy.io.savemat('{}_{}/Unfolding_trajectories.mat'.format(protein, directory), dict(times=times, labels=labels))
    
    
    #TODO: Run Matlab script 
    
    matlab_dict=scipy.io.loadmat('{}_{}/Unfolding_states.mat'.format(protein, directory))
    
    states=matlab_dict['states'][0]-1  #cluster assignmetns for each point
    emit_prob=matlab_dict['emit_prob']
    
    
    
    Nstates=3
    
    clusters=[[] for n in range(Nstates)]
    

    for l, p in enumerate(unique_points):
        assignment=np.argmax(emit_prob[:,l])
        clusters[assignment].append(p)
    
    
    
    
    
    #unique_labels=np.unique(labels)
    #n=len(unique_points)    
    #scipy.io.savemat('Unfolding_trajectories.mat', dict(times=times, labels=labels))
    
    
    
##########################################################################    
    #Here is how I was considering doign distance computation but it didn't seem to work...
    
    
    n=len(unique_points)
    #build transition matrix...am doing this myself rather than using the Markov_transmat function in Folding_rates, since that function is UGLY
    counts=np.zeros((n,n))
    norms=np.zeros((n,n))
    
    for t, i in enumerate(labels):
        if t<len(labels)-1 and times[t+1]!=0:
            i=int(i)
            j=int(labels[t+1])
            counts[i,j]+=1
            norms[i,:]+=1
    transmat=np.divide(counts, norms)
           
    
    #now we compute that distance metric
    
    labels=np.array(labels)
    distance=np.nan*np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            P_i_given_ij=len(np.where(labels==i)[0])/(len(np.where(labels==i)[0])+len(np.where(labels==j)[0]))
            P_j_given_ij=1-P_i_given_ij
            
            distance[i,j]=  P_i_given_ij/transmat[i,j] + P_j_given_ij/transmat[j,i]
            #distance[i,j]=  (transmat[i,j]+transmat[i,i]) *P_i_given_ij/transmat[i,j] + (transmat[j,i]+transmat[j,j]) *P_j_given_ij/transmat[j,i]
            
    distance[np.diag_indices(n)]=0  #sometimes it isn't small...
    
    
    thresh=30
    unused_output, clusters, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(thresh, unique_points, distance, sort_orphans=True, min_clustersize=1)    
   
    
    
    #sort clusters in order of  average nubmer of substrcutreus
    mean_sums=[np.mean([  np.sum( int(l) for l in top ) for top in c]) for c in clusters]
    indices=np.argsort(mean_sums)
    clusters=[clusters[i] for i in indices]
    nclusters=len(clusters)
    unfolding_clusters=np.zeros(len(unfolding_points))
    for p, point in enumerate(unfolding_points):
        unfolding_clusters[p]=[c for c in range(nclusters) if point in clusters[c]][0]
    
    
    
    ###########################################
    #This is the method I was using to do distance computation on July 25, which worked pretty well for full ADK, although maybe less for trunc ADK
    #The idea is that two configurations are defined as kinetically related if there is a comparable number of back and forth exchanges between them

    
    n=len(unique_points)
    
    transitions_count=np.zeros((n,n))
    
    for t in range(len(unfolding_points)-2):
        if unfolding_points[t+1]!=unfolding_points[t] and times[t+1]!=0:
            #if  unfolding_points[t+1]==unfolding_points[t+2]: #we only consider transitions that last more than one timestep
                i=unique_points.index(unfolding_points[t])
                j=unique_points.index(unfolding_points[t+1])
                transitions_count[i,j]+=1
        elif unfolding_points[t+1]==unfolding_points[t] and times[t+1]!=0:
                i=unique_points.index(unfolding_points[t])
                transitions_count[i,i]+=1
    
    distance=np.zeros((n,n))

    
    adjacency=np.zeros((n,n))
    adjacency[np.diag_indices(n)]=1
    for i in range(n):
        for j in range(n):
            pair=[transitions_count[i,j], transitions_count[j,i]]
            if np.min(pair)!=0:
                distance[i,j]=np.max(pair)/np.min(pair)
            else:
                distance[i,j]=10000000
            
    distance[np.diag_indices(n)]=0

    thresh=1.3
    unused_output, clusters, mean_intercluster, mean_intracluster=LoopCluster.loop_cluster_contacts(thresh, unique_points, distance, sort_orphans=True, min_clustersize=1)    
    
    
    #sort clusters in order of  average nubmer of substrcutreus
    mean_sums=[np.mean([  np.sum( int(l) for l in top ) for top in c]) for c in clusters]
    indices=np.argsort(mean_sums)
    clusters=[clusters[i] for i in indices]
    nclusters=len(clusters)
    unfolding_clusters=np.zeros(len(unfolding_points))
    for p, point in enumerate(unfolding_points):
        unfolding_clusters[p]=[c for c in range(nclusters) if point in clusters[c]][0]

    return clusters, unfolding_clusters



#mport sklearn
#from sklearn.hmm import MultinomialHMM

def Histogram(labels, key):
    plt.figure()
    a=plt.hist(labels, bins=len(key)-1)



################################


 
#Beautify the following later  
#Assignments, PDB_files, native_contacts, substructures=joblib.load('FABG_tighter_temps/Substructure_assignments.dat')

#Substructure_scores, PDB_files, native_contacts, substructures=joblib.load('FABG_tighter_temps/Substructure_scores.dat')



#Sort assignments so that they are in order
#PDB_files, N_participants=visualize_PCA.sort_projections(PDB_files, N_participants) 


def Score_snapshot(snapshot, substructures, atom='CA', min_seq_separation=8 ):
    """
    Assigns a set of scores for a snapshot
    the ith score tells you what is the average distnace between pairs of residues residues that participate in the ith substructure, in this snapshto
    If the score is close to the characteristic contact distnace, then the substructure should be mostly formed
    """
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    distances=Analyze_structures.compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
    nsubstructures=np.shape(substructures)[2]
    scores=np.zeros((nsubstructures))
    for s in range(nsubstructures): 
        sub=substructures[:,:,s]
        participation=np.multiply(distances, sub)#gives the overlap between native substrucutres and this snapshot's contacts
        scores[s]=np.mean(participation[np.nonzero(participation)])
    return scores
    

def Count_natives_in_substructures(snapshot, substructures, d_cutoff, atom='CA', formation_cutoff=10, min_seq_separation=8):
    """
    Counts the number of contacts in the substructures for any substructures that
    are "formed", whcih we take to mean that the average distance between the residues
    that participate in that substrucutre is less than formation_cutoff
    """
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    distances=Analyze_structures.compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
    nsubstructures=np.shape(substructures)[2]
    n_natives=0
    n_nonnatives=0
    
    for s in range(nsubstructures): 
        sub=substructures[:,:,s]
        participation=np.multiply(distances, sub)#gives the overlap between native substrucutres and this snapshot's contacts
        score=np.mean(participation[np.nonzero(participation)])
        if score<=formation_cutoff:
            print(s)
            natives_in_sub=np.zeros(np.shape(sub))
            natives_in_sub[np.where((participation<d_cutoff) & (participation!=0))]=1
            n_natives+=np.sum(natives_in_sub)
    return n_natives
        


def Count_native_contacts(snapshot,native_contacts, d_cutoff, atom='CA', min_seq_separation=8, filter_width=2, mode='Hard'):
    """
    Given a snapshot and a matrix of native contacts, counts the number of native contacts in the given
    snapshot, as well as the mean distance between the residue pairs that participate in native contacts
    
    d_cutoff is the maximum allowed distnace between two residues for them to be considered a contact when mode is 'Hard'
    If mode is 'Soft', applies a sigmoidal function to define contacts
    """
    coords, resis=Analyze_structures.read_PDB(snapshot, atom)
    distances=Analyze_structures.compute_contacts_matrix(coords, mode='distances', min_seq_separation=min_seq_separation)
    filtered_distances=np.multiply(distances, native_contacts)
    mean_distance=np.mean(filtered_distances[np.nonzero(filtered_distances)])
    natives_in_snapshot=np.zeros(np.shape(filtered_distances))
    if mode=='Hard':
        natives_in_snapshot[np.where((filtered_distances<d_cutoff) & (filtered_distances!=0))]=1
        n_natives=np.sum(natives_in_snapshot)
    else:
        natives_in_snapshot[np.where(filtered_distances!=0)]=d_cutoff/(d_cutoff+filtered_distances[filtered_distances!=0])
        n_natives=np.sum(natives_in_snapshot)
    
    
    #We want to compute non-native contacts without caring too much about native contacts that have been sligthly shifted
    #To do this will eliminate all contacts that are within filter_width of the native contacts
    

    native_filter=cp.deepcopy(native_contacts)
  
    for i in range(-filter_width, filter_width+1):
        native_filter+=np.roll(native_contacts, i, axis=0)
        native_filter+=np.roll(native_contacts, i, axis=1)

    native_filter[np.where(native_filter)]=1
    
    
    nonnative_distances=np.multiply(distances, 1-native_filter)  #distances between residues that don't interact in native structure
    nonnatives_in_snapshot=np.zeros(np.shape(nonnative_distances))
    nonnatives_in_snapshot[np.where((nonnative_distances<d_cutoff) & (nonnative_distances!=0))]=1
    n_nonnatives=np.sum(nonnatives_in_snapshot)
    
    return n_natives, n_nonnatives
        



def Substructure_scores_to_barcodes(data, gmix=True, thresh=10 ):
    """
    Input an array of substructure scores, data
    By default, uses a Gaussian mixture model to learn two clusters
    """
    barcodes=[]
    if gmix:
        barcodes_array=np.zeros(np.shape(data))  
        for i in range(np.shape(data)[1]):
            x=data[:,i]
            indices=np.where(x<30)
            x=x[indices]  #eliminate long tails that may confuse the algorithm
            gmix, labels, ll = ClusterPCA.Learn_clusters(x, 2, n_trials=1, training_set=None, dim=1)
            gmix, barcodes_array[indices,i] = ClusterPCA.fix_cluster_order(gmix, labels, ClusterPCA.Determine_autosort_order(gmix, reverse=True))
        for j, y in enumerate(barcodes_array):
            string=''
            for s in y:
                string='{}{}'.format(string, int(s))
            barcodes.append(string)
    else:
    
        for i in range(len(data)):
            string=''
            for j in range(len(data[i,:])):
                if data[i,j]<=thresh:
                    string='{}1'.format(string)
                else:
                    string='{}0'.format(string)
            barcodes.append(string)
    #labels, key = Bernie_elimination(barcodes, PDB_files, 0, plot=False)
    return barcodes
    


def Bernie_elimination(Assignments, PDB_files, Sanders_thresh, plot=True):
    """
    Eliminate assignments that do not represent more than some fraction of the population given by Sanders_thresh
    """    
    
    #Assign a numerical label to each assignment in roughly order of progressive unfolding
    #
    unique=list(set(Assignments))
    
    #We whipe out clusters that represent too small a fraction of the populaiotn
    
    
    
    Bernie_dic={} #this dic will tell you what fraction of the population is represented by a given cluster 
    for u in unique:
        Bernie_dic[u]=len([i for i in Assignments if i==u])/len(Assignments)
    
    
    
    for i, assign in enumerate(Assignments):
        if Bernie_dic[assign]<Sanders_thresh:
            Assignments[i]='BustBigBank' #label these for removal
    
    
    unique=list(set(Assignments))
    if 'BustBigBank' in unique:
        unique.remove('BustBigBank')
    
    sums=[-sum([int(i) for i in str]) for str in unique]
    sort=natsort.natsorted(zip(sums, unique))
    dic={}  #a dictionary to readily map assignments to labels 
    
    sum_dic={} #total numbe of substructures 
    
    key=[]
    for i in range(len(sort)):
        key.append(sort[i][1])
        dic[sort[i][1]]=i
        
        sum_dic[sort[i][1]]=sort[i][0]
    dic['BustBigBank']=np.nan
        
    labels=[dic[a] for a in Assignments] 
    if plot:
        Histogram(labels, key)

    return labels, key


#as an alternative, we assign label 1 if and only iff the amount of substructures is less than or equal to thresh
#thresh=3
#labels=[]


#for a in Assignments:
#    if sum_dic[a]>=-thresh:
#        labels.append(1)
#    else:
#        labels.append(0)




