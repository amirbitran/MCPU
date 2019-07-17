# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:11:38 2017


Filters PCA data and clusters it
@author: amirbitran
"""

import joblib
import numpy as np
#import matplotlib.pyplot as plt
import natsort  #used to sort things
#from mpl_toolkits.mplot3d import Axes3D
import sklearn
from sklearn import cluster
import fnmatch
import visualize_PCA
from sklearn import mixture
import pickle
import networkx as nx
import scipy
from scipy import spatial
from hmmlearn import hmm
import copy as cp
import Find_identical_trajectories
import ClusterSubstructures

def dist(p1, p2):
    return np.sqrt(np.dot(p1-p2, p1-p2))



def pairwise_distance_matrix(coords):
    """
    Constructs a pairwise distance matrix for all points in np array coords
    """
    p_dist=np.zeros((len(coords),len(coords)))
    for i in range(len(coords)):
        for j in range(i):
            p_dist[i,j]=dist(coords[i,:], coords[j,:])
    return p_dist+np.transpose(p_dist)



def Filter_points(projections, PDB_files, trajectories, epsilon=8, required_neighbors=2):
    """
    Separates out unfolded snapshots using the observation that unfolded trajectories tend to be highly random in PCA space
    For each trajectory, we insist that every point lie within distance epsilon of at least some number of other points in that trajectory
    that number is given by required_neighbors
    If this condition is satisfied, then the point is kept, as we assme it is sitting at a free energy minimm
    Otherwise, it is discarded
    
    This function does not seem to work super well from what I remember
    
    """
    
    folded_coords=[]
    folded_files=[]
    
    unfolded_coords=[]
    unfolded_files=[]
    
    #First, we loop through the trajectories we care about
    for curr_traj in trajectories:
        files, coords = visualize_PCA.Get_trajectory(projections, PDB_files,visualize_PCA.num2str(curr_traj))
        
        #Next, we construct a pairwise distance matrix for all points
        p_dist=pairwise_distance_matrix(coords)
        
        #Now, for each point p, we count the number of points whose distance from p is less than epsilon
        
        for n,p in enumerate(coords):
            n_neighbors=len(np.where(p_dist[n,:]<epsilon)[0])-1   #we subtract 1 because trivially, the point will be a distance zero from itself
            if n_neighbors>=required_neighbors:
                folded_coords.append(p)
                folded_files.append(PDB_files[n])
            else:
                unfolded_coords.append(p)
                unfolded_files.append(PDB_files[n])
    visualize_PCA.ScatterPlot2D(np.vstack(folded_coords))
    return folded_coords, folded_files, unfolded_coords, unfolded_files
                
                


def Learn_clusters(projections, nclusters, n_trials=10, training_set=None, dim=2):
    """
    Runs Gaussian mixture model multiple times and returns parameter set that fits all data with highest log likelihood
    You have the option of learning parameters using only a subset of data specificed by training_set.    
    I know this can supposedly be done directly using sklearn by setting n_init, but this doesn't seem to work for me...changing n_init (even up to 100) still yields local maxima
    """
    ll=[]
    trials=[]
    projections=np.atleast_2d(projections)
    
    if np.shape(projections)[0]==1:
        projections=np.transpose(projections)
    for trial in range(n_trials):
        if training_set!=None:
            traj_files, traj_coords = visualize_PCA.Get_trajectory(projections, PDB_files,training_set)
            gmix=sklearn.mixture.GaussianMixture(n_components=nclusters, n_init=n_trials).fit(traj_coords[:,0:dim])
        else:
            gmix=sklearn.mixture.GaussianMixture(n_components=nclusters, n_init=n_trials).fit(projections[:,0:dim])
        
        ll.append(gmix.score(projections[:,0:dim]))
        trials.append(gmix)
    best_ll=max(ll)
    best_trialnum=[l for l in range(len(ll)) if ll[l]==best_ll][0]
    best_gmix=trials[best_trialnum]
    labels=best_gmix.predict(projections[:,0:dim])
    return best_gmix, labels, best_ll
    
def ScatterClusters(projections, labels, ax1title='rmsd', ax2title='contacts'):
    plt.figure()
    if type(labels)==list: labels=np.array(labels)
    cmap = plt.get_cmap('Paired')
    nclusters=len(set(labels))
    colors = [cmap(i) for i in np.linspace(0, 1, nclusters)]
    for c in range(nclusters):
        #locs=np.where(kmeans.labels_==c)[0]
        locs=np.where(labels==c)[0]
        plt.scatter(projections[locs, 0], projections[locs,1], color=colors[c], label='Cluster {}'.format(c))
        #ax.scatter(projections[locs, 0], projections[locs,1],projections[locs,2], color=colors[c])
   
    
    plt.xlabel(ax1title, fontsize=25)
    plt.ylabel(ax2title, fontsize=25)
    plt.legend()
    plt.title('Gaussian mixture clustering', fontsize=25)  
    plt.tick_params(axis='both', which='major', labelsize=22, pad=2)
    
def Cluster2D(projections, nclusters, PDB_files, training_set=None):
    """
    Clusters PCA projections into nclusters using EM mixture model
    Also have the optin of clustering only the files in a given set of trajectories given by training_set...
    these files are then used as the training set to construct the clusters
    
    Only first two dimsnesions of data are used
    """
    gmix, labels, best_ll=Learn_clusters(projections, nclusters, training_set=training_set)
    print('The optimal log likelihood of all data given model is {}'.format(best_ll))
    #fig=plt.figure()
    #ax=fig.add_subplot(111, projection='3d')
    #colors=['r','b','k','m', 'g', 'y', 'c']

    plt.figure()
    for c in range(nclusters):
        #locs=np.where(kmeans.labels_==c)[0]
        locs=np.where(labels==c)[0]
        plt.scatter(projections[locs, 0], projections[locs,1], label='Cluster {}'.format(c))
        #ax.scatter(projections[locs, 0], projections[locs,1],projections[locs,2], color=colors[c])
    
    plt.xlabel('Principal component 1')
    plt.ylabel('Principal component 2')
    plt.title('Gaussian mixture clustering')  
    plt.legend()
    return gmix, labels






def Hidden_markov(gmix, projections, PDB_files, prior='Rare', diag=100, upper_triang=10):
    """
    prior can take on any of the following values:
    
    Rare: We make all transitions very unliekly
    Irreversible: We make transitions from higher to lower clusters very unliekly
        In this case, we initialize transmat dirichlet with diagonal entries of value diag
        and upper triangular (excluding diagonal) entries of value upper_triang
     None: No prior
    """
    
    
    traj_len=200
    means=gmix.means_
    
    if type(gmix)==hmm.GaussianHMM:
        covariances=gmix.covars_
    else:
        covariances=gmix.covariances_
    NStates=np.shape(covariances)[0]
    
    files, trajectories=visualize_PCA.Get_trajectory(projections, PDB_files, 'unfolding')
    lengths=[traj_len for x in range(int(np.shape(trajectories)[0]/traj_len))]
    
   
         
    startprob=np.zeros(NStates)
    startprob[0]=1  #everybody starts folded
    
    if prior=='Rare':
         transmat_prior=np.ones((NStates, NStates))
         transmat_prior[np.diag_indices(NStates)]=10000000000
         #transmat_prior[np.diag_indices(NStates)]=10
         model=hmm.GaussianHMM(n_components=NStates, covariance_type="full", transmat_prior=transmat_prior, params='t', init_params='t')
    elif prior=='Irreversible':  
         transmat_prior=np.ones((NStates, NStates))
         #transmat_prior[np.tril_indices(NStates)]=1
         transmat_prior[np.triu_indices(NStates)]=upper_triang
         transmat_prior[np.diag_indices(NStates)]=diag
         
         model=hmm.GaussianHMM(n_components=NStates, covariance_type="full", params='t', init_params='t', transmat_prior=transmat_prior)
    elif prior==None:
        model=hmm.GaussianHMM(n_components=NStates, covariance_type="full", params='t', init_params='t')
    else:
        print('Prior can either be "Rare", "Irreversible", or None !!!')
    model.startprob_ = startprob
    model.means_ = means
    model.covars_ = covariances
            
    model.fit(trajectories, lengths)

    return model

    
def Determine_autosort_order(gmix, reverse=False):
    #M=gmix.means_[:,0:-1]  #exclude natives, which if exists is usually at the end
    M=gmix.means_
    if reverse: M=-M
    mean_magnitudes=[]
    for clust in M:
        mean_magnitudes.append(np.sum(clust))  #Most reaction coordinates tend to increase as protein becomes more unfolded
        #Thus we sort in order of increasing sum of reaction coordinates
        
        #Old method:
        #mean_magnitudes.append( np.sum(clust[np.where(clust>0)]**2)) #only compute magnitudes for positive values
        #We want to sort by magnitude, but if we include negative values, this gets confused since sometimes the most folded clusters have very negative values
    return np.argsort(mean_magnitudes) 
    


def fix_cluster_order(gmix, labels, order):
    """
    Suppose the ith entry of order has value j
    Then this tells you that entries which should have been assigned to label i have currently been assigned cluster j
    Note that this is opposite of how swap_labels works

    In swap labels:    
    Suppose entry i of order has value j. Then each time label i appears in the original assignments, it is replaced with label j

    """
    M=gmix.means_
    C=gmix.covariances_
    NStates=np.shape(C)[0]
    new_labels=cp.deepcopy(labels)
    
    
    means=np.zeros(np.shape(M))
    covariances=np.zeros(np.shape(C))
    
    for i in range(NStates):
        means[i,:]=M[order[i],:]
        covariances[i,:,:]= C[order[i],:,:]
        new_labels[np.where(labels==order[i])]=i
        
         
    gmix.means_ = means
    gmix.covariances_ = covariances
    
    return gmix, new_labels
    
    

def Split_cluster(cluster, data, labels, PDB_files, nclusters):
    """
    THIS DOES NOT WORK!!!
    
    Takes all snapshots assigned to a given cluster (ex. 0,1,2) and applies a Guassian
    mixture algorithm just to these snapshots, thus splitting an existing cluster into further clusters
    """
    #First, get all PDB files assigned to cluster
    indices_of_interest=[]
    data_of_interest=[]
    
    for n,f in enumerate(PDB_files):
        if labels[n]==cluster:
            indices_of_interest.append(n)
            data_of_interest.append(data[n,:])
            
    data_of_interest=np.array(data_of_interest)

    gmix, new_labels, best_ll=Learn_clusters(data_of_interest, nclusters, training_set=None, dim=np.shape(data)[1], n_trials=10)
    
    
    #We assign labels as follows:
    #First, order the new lables in order of ascending means, as usual
    #Then, the snapshots wiht the highest label are given the original label value "cluster"
    #Remaining snapshots are kept in order, but their values are shifted down so that they do not conlfict with exisiting lables
    order=Determine_autosort_order(gmix)
    gmix, new_labels=fix_cluster_order(gmix, new_labels, order)
    
    
    maxlabel=np.max(new_labels)
    shift=(maxlabel-1)-np.min(labels)+1
    for i, l in enumerate(new_labels):
        n=indices_of_interest[i]
        if l!=maxlabel:
            labels[n]=labels[n]-shift
            
    
    
    
    
    
    



def Cluster(projections, nclusters, PDB_files, training_set=None, HMM=True, HMM_mode='Multi', autosort=True, prior='Irreversible'):
    """
    Clusters PCA projections into nclusters using EM Gaussian mixture model
    Also have the optin of clustering only the files in a given set of trajectories given by training_set...
    these files are then used as the training set to construct the clusters
    
    All dimensions of data are used
    If autosort=True, then function automatically orders cluster labels based on the first attribute of the cluster means, in increasing order
    
    You have the option of using or not using a HMM
    If you do use a HMM, there are two modes:
    If HMM_mode='Single', you fit all the unfolding trajectories. You then return a single transition matrix
    Otherwise if HMM_mode='Multi', you do a HMM fit on each individual temperature. In that case you return a set of transition matrix stored in transmat, each page of which is the transition matrix at a given temp

    """
    #projections=np.atleast_2d(projections)    
    gmix, labels, best_ll=Learn_clusters(projections, nclusters, training_set=training_set, dim=np.shape(projections)[1], n_trials=10)
    print('The optimal log likelihood of all data given model is {}'.format(best_ll))
    
    
    #Re-sort ordering of clusters
    if autosort==True:
        order=Determine_autosort_order(gmix)
    else:
        print("Need to fix function so that custom orderings can be inputted!" )
    gmix, labels=fix_cluster_order(gmix, labels, order)

    if HMM==True:
        print('Training hidden Markov model on unfolding simulations' )
        
        if HMM_mode=='Single': #we make this a simple HIdden markov model trained at all temps
            gmix=Hidden_markov(gmix, projections, PDB_files, prior=prior)
        elif HMM_mode=='Multi':
            #First, obtain a list of all unique temperatures at which unfolding simulations were run 
            unfolding_files, unfolding_trajectories=visualize_PCA.Get_trajectory(projections, PDB_files, '*/*/*.***_*.*' ) #all unfolding files
            temps=[f.split('_')[-2] for f in unfolding_files]
            temps=list(set(temps))
            
            #sort the temperatures
            sort=natsort.natsorted(zip(temps, [float(T) for T in temps]))
            temps=[tup[0] for tup in sort]
            
            n_temps=len(temps)
            
            transmat=np.zeros((nclusters, nclusters, n_temps))
            
            #Consider doing the following:
            gmix.temps=temps              
            
            for n, temp in enumerate(temps):
                f, trajectories = visualize_PCA.Get_trajectory(projections, PDB_files, '{}_'.format(temp) )
                hm=Hidden_markov(gmix, trajectories, f, prior=prior)
                transmat[:,:,n]=hm.transmat_

        #I think the following is perfectly legal to do! That way you can return the transmat as part of gmix, instead of as a separate variale
            gmix.transmat_=transmat 
    
        
    if np.shape(projections) [1]>1:    
        ScatterClusters(projections, labels)    
    return gmix, labels

def lookup_labels(file, PDB_files, labels):
    """
    returns the label for a given file or set of files
    file is a string that is contained within a single file, or a set of files
    (ex. '720.199000', or '0.850_15')
    """
    indices=[f for f in range(len(PDB_files)) if fnmatch.fnmatch( PDB_files[f], '*{}*'.format(file)) ]
    for index in indices:
        print("{}: {} \n".format(PDB_files[index], labels[index]))
    #return [labels[index] for index in indices]
     
def lookup_files(label,labels, PDB_files, traj='All'):
    """
    returns all files that are assigned to a given label   
    You can use traj to specify that you only care about files in a given trajectory
    The defalt is 'All' (all trajectories are considered )
    """
    files=[f for n,f in enumerate(PDB_files) if labels[n]==label ]
    if traj!='All':
         files = [f for f in files if fnmatch.fnmatch( f, '*_{}*'.format(traj))]
    return files
         
def Elbow_test(projections, training_set=None,  cluster_range=np.arange(1,8,1)):
    """
    runs clustering using various numbers of clusters, and determines log likelihood that the clusters produce the FULL data
    As before, you have the option of training the algorithm on a partial dataset given by training_set    
    """
    ll=[]
    heritages=[]  #keeps track of what fraction of files in some cluster for n_clusters belonged to some other cluster for n_clusters-1
    for nclusters in cluster_range:
        if nclusters>1: old_labels=labels
        print('Computing optimal LL for {} clusters'.format(nclusters))
        gmix, labels, curr_ll=Learn_clusters(projections, nclusters, training_set=training_set, dim=np.shape(projections)[1], n_trials=10)  
        order=Determine_autosort_order(gmix)
        gmix, labels=fix_cluster_order(gmix, labels, order)
        ll.append(curr_ll)
        
        
        #now, we build a matrix called heritage that traces the heritage of the files in the current cluster. This is an n_clusters-1 by n_clusters
        #matrix in which entry (i,j) tells you, for those files in cluster j when we have n_clusters, what fraction came from cluster i when we had
        #n_clusters-1? Thus, columns sum to 1
        
        if nclusters>1:
            heritage=np.zeros((nclusters-1, nclusters))
            for j in set(labels):
                files_with_label_j=np.where(labels==j)
                parents=old_labels[files_with_label_j]
                for i in set(old_labels):
                    files_j_descendedfrom_i=np.where(parents==i)[0]
                    heritage[i,j]=len(files_j_descendedfrom_i)/len(files_with_label_j[0]) #tells you what fraction of snapshots that currently belong to cluster j were previously in cluster i
            heritages.append(heritage)
            
    plt.figure()
    plt.plot(cluster_range,ll)
    plt.title('Log likelihood')
    plt.xlabel('Number of clusters')
    
    Visualize_lineage(heritages)
    return heritages
        

def Predict_unfolding_clusters(sim_dir, model_dir, traj, mode='MLE', combine_labels=None, filter_osc=True, apply_correction=False, plot=True):    
    """
    THIS IS OUT OF DATE!!
    Predicts optimal clusters for unfolding simulation trajectories
    Loads unfolding simulation data stored in sim_dir and obtains trajectory traj
    Then loads Gaussian mixture model stored in model_dir
    Finally, uses parameters of model_dir to find best fit for trajectory traj
    
    mode tells you how you assign snapshots to cluster. If mode='distance', you assign each PC projection to the nearest centroid. If mode='MLE', assign to the cluster most likely to produce that point ina Gaussian mixture model
    
    ex. Predict_unfolding_clusters('ADK_unfolding/PC_projections_T_0.800.dat', 'ADK_tighter_temps', '0.800_3')
    predicts the clusters for unfolding simulation whose projections onto PCs are stored in ADK_unfolding/PC_projections_T_0.800.dat
    using the model parameters stored in ADK_tighter temps, for trajectory 0.800_3
    
    NOTE: Changed this so that, for model_dir, you also need to specify a .dat file
    
    Note also that, for some proteins, you may have artifically combined labels into one, ex. 0 and 1.
    The problem is, this does not change the internal parameers of the Gaussian mixture model.
    So the Gaussian mixture model makes its predictions with the old labelign system, but you need to re-specify which labels you want to combine post-hoc
    """

    [ pca, projections, PDB_files]=visualize_PCA.load_PCA(sim_dir)
    [ traj_files, traj_coords]=visualize_PCA.Get_trajectory(projections, PDB_files, traj)
    
    [gmix, labels, PDB_files] = joblib.load(model_dir)
    if mode=='distance':
        centroids=gmix.means_   #n_clusters by n_dimensions
        distances=scipy.spatial.distance.cdist(traj_coords, centroids)  #index i,j of this matrix tells you the distances between trajectory projection i and centroid j
        sim_labels=np.argmin(distances, axis=1)
    elif mode=='MLE':
        sim_labels=gmix.predict(traj_coords[:,0:len(gmix.means_[0,:])])
    else: 
        print('mode must be either "distance" or "MLE"')
    sim_probs=gmix.predict_proba(traj_coords[:,0:len(gmix.means_[0,:])])  #posterior probabilities that each datapoint (inexed by rows) is produced by each cluster (indexed by columns)
    likelihoods=Gaussian_scores(traj_coords, gmix)
    total_scores=gmix.score_samples(traj_coords[:,0:len(gmix.means_[0,:])])
    if filter_osc==True: sim_labels=filter_oscillations(sim_labels)
    if apply_correction==True: 
        problem_cluster=2
        suspect_correct_cluster=1
        sim_labels=Correct_clustering(sim_labels, sim_probs, problem_cluster, suspect_correct_cluster)
    #print(total_scores)
    if combine_labels!=None:
        sim_labels=Combine_clusters(combine_labels[0], combine_labels[1], sim_labels)
    
    if plot==True:
        plt.figure()
        plt.plot(sim_labels)
        plt.ylim(0, np.max(labels)+0.25)
        
        
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        plt.title('Log posterior probabilities of clusters given data')
        plt.xlabel('MC step (x 1,000,000)')
        for n,col in enumerate(sim_probs.transpose()):
            ax.plot(col, label='Cluster {}'.format(n))
        ax.set_yscale('log')
        ax.set_ylim([10**(-15), 2])
        plt.legend(loc='lower left')
    

    
    return sim_labels, traj_files, likelihoods



def filter_nans(sim_labels, mode='prev'):
    """
    Wherever you see a stream of Nan's, do one of two things:
    
    1. if mode='prev', replace everything in that stream with the last non-Nan value that came before
    2. If mode='next', replace everythign in that stream with the first non-Nan value that comes next
    """
    for n,l in enumerate(sim_labels):
        if np.isnan(l):
            if mode=='prev':
                if n>0:  #you have non-nan timesteps before this to work with
                    sim_labels[n]=sim_labels[n-1]
                else: #need to switch to "next" mode
                    remaining=sim_labels[n:]
                    sim_labels[n]=remaining[np.min(np.where(~np.isnan(remaining))[0])]
                    
                
                
                
            elif mode=='next': 
                remaining=sim_labels[n:]
                if len(np.where(~np.isnan(remaining))[0])==0: #In this case, everything from here to the end is Nan, so we go with what came before
                    sim_labels[n:]=sim_labels[n-1]
                else:
                    sim_labels[n]=remaining[np.min(np.where(~np.isnan(remaining))[0])]
    return sim_labels
            
        


def filter_oscillations(sim_labels, thresh=7, islands_only=True):
    """
    Gets rid of oscillations in unfolding trajectory that are shorter than thresh
    That is, this code considers every "island" of snapshots comprised of one label, where the surrounding snapshots
    are one other label.
    If this island is smaller than thresh, then the labels for snapshots within this island are changed
    to correspond to the surrounding labels
    
    By default, we only filter out islands.
    But if islands_only=False, we also filter out rapid changes between indices, ex. 
     000112222. In this case  the code gets rid of the 
    middle index (1) so long as the last index (2) hangs around for longer than thresh
    """
    changes=[]    #Keeps track of indices in simulation at which change occurs
    #accepted_changes=[]
    change_n=-1  #Enumerates the changes that have occured
    #accepted_change_n=-1
    
    
    sim_labels=filter_nans(sim_labels) #first, get rid of any nans that may be present
    
    if islands_only==True:
        for n in range(len(sim_labels)):  #loop through all points in simulation
            if n>0 and sim_labels[n]!=sim_labels[n-1]:  #There is a change between time n-1 and time n
                changes.append(n)
                change_n+=1
                if change_n>0 and changes[change_n]-changes[change_n-1]<thresh and sim_labels[changes[change_n]]==sim_labels[changes[change_n-1]-1]:   #We spent too little time between the last two changes, so we  filter what comes in between the changes
                    sim_labels[changes[change_n-1]:changes[change_n]]=sim_labels[n]  #Filter out what comes between the changes
                    changes.remove(changes[change_n])  #In this case, we remove the last two changes, since we actually stay in the same place in the period between these 'changes'
                    changes.remove(changes[change_n-1])
                    change_n=change_n-2  #we eliminated the last two changes, so reduce our change counter by 2

                        
            if len(changes)>0:
                if n==len(sim_labels)-1 and n - changes[-1]<thresh:  #If, after all this filtering, the last few labels are different than what came before (specificially, if we have less than thresh of these differnet labels at the end), then we do not accept that change in label since we do not have enouhg points to tell whehter is a fluctuation or real
                    sim_labels[changes[-1]:len(sim_labels)]=sim_labels[changes[-1]-1]
    else:
        for n in range(len(sim_labels)):  #loop through all points in simulation
            if n>0 and sim_labels[n]!=sim_labels[n-1]:  #There is a change between time n-1 and time n
                changes.append(n)
                change_n+=1
                if change_n>0 and changes[change_n]-changes[change_n-1]<thresh:   #We spent too little time between the last two changes, so we  filter what comes in between the changes
                    sim_labels[changes[change_n-1]:changes[change_n]]=sim_labels[n]  #Filter out what comes between the changes
                    if sim_labels[changes[change_n]]==sim_labels[changes[change_n-1]-1]:  #We filtered out an "island": a set of indices predicted to be on cluster while the surrounding indices are all one other cluster (ex. 22222 1 22222, we replace the lone 1 with a 2)
                        changes.remove(changes[change_n])  #In this case, we remove the last two changes, since we actually stay in the same place in the period between these 'changes'
                        changes.remove(changes[change_n-1])
                        change_n=change_n-2  #we eliminated the last two changes, so reduce our change counter by 2
                    else:  #We filtered out a change that is not an island--for instance, if 0 goes to 2 very briefly then to 1. 
                        changes.remove(changes[change_n])  #In this case, we remove the current change, since we have rejected it--we just care about the previous change position
                        change_n=change_n-1
                        
            if len(changes)>0:
                if n==len(sim_labels)-1 and n - changes[-1]<thresh:  #If, after all this filtering, the last few labels are different than what came before (specificially, if we have less than thresh of these differnet labels at the end), then we do not accept that change in label since we do not have enouhg points to tell whehter is a fluctuation or real
                    sim_labels[changes[-1]:len(sim_labels)]=sim_labels[changes[-1]-1]
             
    return sim_labels
                    

def score_trajectory(directory, traj):
    """
    You input a directory and a trajectory, this returns a list of scores at each time point.
    The score at a given time point is the total log likelihood that the time point is produced by any of the clusters in
    the Gaussian mixture model stored in that directory
    """
    [ pca, projections, PDB_files]=visualize_PCA.load_PCA('{}/PCA_complete.dat'.format(directory))
    [traj_files, traj_coords]=visualize_PCA.Get_trajectory(projections, PDB_files, traj)
    [gmix, labels, PDB_files] = joblib.load("{}/Gaussian_labels.dat".format(directory))
    total_scores=gmix.score_samples(traj_coords[:,0:len(gmix.means_[0,:])])
    return total_scores
    


def Correct_clustering(sim_labels, sim_probs, problem_cluster, suspect_correct_cluster, thresh=-60):
    """
    Inspired by the observation that, for full ADK unfolding trajectories, clustering algorithm tends to incorrectly assign snapshots to an
    partially unfolded cluster (cluster 2), when in fact, the snapshots did not have the characteristic of that partially unfolded cluster and 
    should be assigne to a mostly folded cluster (cluster 1)
    I also observed that these incorrect assignments tend to occur when the posterior probability of cluster 1 hovered around 10^-60, sometimes spiking up
    to much higerh values. WHen that posterior probability was less than 10^-100 (and sometimes 0, below the precision limit), the assignment to cluster 2 tended to be correct to my eye
    
    This function corrects these sorts of mistakes in a very heuristic way. You input a series of assignments (given by lables), which have generally been filtered
    using filter_oscillations, as well as your posterior probabilities for all the clusters at each time point. You also input a problem cluster label that tends to be
    incorrectly asigned (ex. problem_clsuter = 2), and the cluster number that you suspect is correct (ex. suspect_correct_cluster=1).
    
    The function then finds all points that have been assigned to problem_cluster and computes the mean log likelihood that those snapshots belong to suspect_correct_cluster
    If that mean log likelihood is greater than thresh, then all thos eassignments are changed to suspect_correct_cluster
    Otherwise, problem_cluster is kept as the assignment
    """
    
    problem_labels=np.where(sim_labels==problem_cluster)
    probs_of_suspect=sim_probs[problem_labels, suspect_correct_cluster]
    log_probs_of_suspect=np.log10(probs_of_suspect[np.where(probs_of_suspect)])  #don't want 0 log probabilities
    mean_log_prob=np.mean(log_probs_of_suspect)
    print(mean_log_prob)
    
    if mean_log_prob>thresh: sim_labels[problem_labels]=suspect_correct_cluster
    return sim_labels
    

def Visualize_lineage(heritages):
    """
    Lets you visualize a matix heritages produced by elbow test
    """
    fig=plt.figure()
    ax=fig.add_subplot(111)
    G=nx.Graph()
    y_coords=np.arange(1,len(heritages)+2,1)
    plt.ylim(0.5, len(heritages)+1.5)
    plt.ylabel('Number of clusters')
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')
    plt.gca().invert_yaxis()
    
    
    pos={}
    for N in range(len(heritages)+1):
        y=y_coords[N]
        x_coords=np.linspace(-N,N,N+1)

        for j,x in enumerate(x_coords):
            G.add_node("N={} \n cluster {}".format(N+1,j), pos=(x,y))
            nx.draw_networkx_nodes(G, nodelist=["N={} \n cluster {}".format(N+1,j)], with_labels=True, node_size=1000, pos={"N={} \n cluster {}".format(N+1,j):(x,y)})
            pos["N={} \n cluster {}".format(N+1,j)]=(x,y)        
            #positions.append((x,y))
            if N>0:
                for i in range(N):
                    W=heritages[N-1][i,j]
                    #print(N,(i,j),W)
                    #print(i)
                    G.add_edge("N={} \n cluster {}".format(N,i), "N={} \n cluster {}".format(N+1,j), weight=W)
                    nx.draw_networkx_edges(G,    pos,      edgelist=[("N={} \n cluster {}".format(N,i), "N={} \n cluster {}".format(N+1,j))], width=W, edge_color=[100*W] ) 
                    #fig.colorbar(ax)
                    
                    #if 0<W<1.0: nx.draw_networkx_edge_labels(G,    pos,      edge_labels={("N={} \n cluster {}".format(N,i), "N={} \n cluster {}".format(N+1,j)) : '{}'.format(np.round(W, decimals=2)) } ) 
                    
                    #print("N={} cluster {} to N={} cluster {} has weight {}".format(N,i,N+1,j,W))


def swap_labels(labels, new_assignments):
    """
    Lets you reassign label numberings
    Input the original lables, as well as new_assignments, which is a list whose length must be equal to the nubmer of distinct original labels
    new_assignments must be some permutation of the unique distinct original labels
    Suppose entry i of new_assignments has value j. Then each time label i appears in the original assignments, it is replaced with label j
    
    FOr instance, suppose you have 4 clusters. You are happy with assignments 0 and 3, but want to swap all assignments of 1 for 2, and vice versa.
    Then you input new_assignments=[0, 2, 1, 3]
    """
    new_labels=[]
    for l in labels:
        new_labels.append(new_assignments[l])
    return new_labels            

def Gaussian_scores(traj_coords, gmix):
    """
    Evaluates the Gaussians specified by gmix at each point in the matrix traj_coords to yield the likelihoods at these points
    Returns scores, which is a matrix whose rows are points whose probabilities we care about, and whose columns are gaussians
    """        
    scores=np.zeros((np.shape(traj_coords)[0], np.shape(gmix.means_)[0]))
    for i in range(np.shape(gmix.means_)[0]):
        scores[:,i]=scipy.stats.multivariate_normal.pdf(traj_coords, gmix.means_[i,:], gmix.covariances_[i,:,:])
    return scores


def Markov_prediction(trajectories, files, gmix, filter_mode='First commit'):
    """
    Obtains unfolding label predictions using hidden markov model gmix,
    then filters by enforcing that all transitions be irreversible
    using some method described in enforce_irreversibiltiy
    
    filter_mode can also be None, in which case you do not enforce irreversibiltiy
    """
    traj_len=np.shape(trajectories)[0]
    #files, trajectories=visualize_PCA.Get_trajectory(data, PDB_files, traj)
    lengths=[traj_len for x in range(int(np.shape(trajectories)[0]/traj_len))]
    
    
    temp=files[0].split('_')[-2] #what temperature is this trajectory at?
    
    if type(gmix)==sklearn.mixture.gaussian_mixture.GaussianMixture: #need to build into a Gaussian mixure model
        n=[i for i in range(len(gmix.temps)) if gmix.temps[i]==temp][0]
        transmat=gmix.transmat_[:,:,n]
        
            
        NStates=np.shape(gmix.covariances_)[0]
        
        startprob=np.zeros(NStates)
        startprob[0]=1  #everybody starts folded
             
        model=hmm.GaussianHMM(n_components=NStates, covariance_type="full")
        
        model.startprob_ = startprob
        model.means_ = gmix.means_
        model.covars_ = gmix.covariances_
        model.transmat_=transmat
        labels=model.predict(trajectories, lengths)
        
    else: #already a Gaussian mixtrue model
        labels=gmix.predict(trajectories, lengths)
    
    if filter_mode!=None: 
        labels=enforce_irreversibility(labels, filter_mode=filter_mode)
        labels=labels[0,:]
    return labels



def enforce_irreversibility(labels, filter_mode='First commit'):
    """
    Forces trajectories to be irreversible in one of two ways
    Either:
    1. You cannot go back to labels you have already visitied (filter_mode='First commit')
    
    or 
    
    2. You do not accept a new label until you stop returning to labels you had visited before that
    (filter_mode='Last commit')
    So for instance, [0,1,2,1,2] would become [0,1,1,1,2], whereas [0,1,2,1,1] would become [0,1,1,1,1]
    
    """
    labels=np.atleast_2d(labels)
    #change_indices=[]  #where do changes occur?
    for traj in labels:
        prev_x=traj[0]
        visited_already=[prev_x]
        for n in range(1, len(traj)):
            x=traj[n]
            if x!=prev_x:
                #change_indices.append(n)
                if x in visited_already:
                    if filter_mode=='First commit':
                        traj[n]=prev_x
                        x=prev_x
                    elif filter_mode=='Last commit':
                        #FIND INDEX FOR LAST TIME YOU WERE AT X...
                        last_visit=np.max(np.where(traj[0:n]==x))
                        #REMOVE ALL THOSE VALUES FROM VISITED_LAREADY
                        for y in set(traj[last_visit+1:n]):
                            visited_already.remove(y)
                        #visited_already.remove([y for y in set(traj[last_visit+1:n])])
                        #REVERT EVERYTHING BETWEEN THAT LAST INDEX AND THE CURENT INDEX SO THAT IT EQUALS X
                        traj[last_visit:n]=x

                else:
                    visited_already.append(x)
            prev_x=x
    return labels
          

    
    

def Plot_Markov_trajectory(gmix, data, PDB_files, traj, filter_mode='First commit'):
    '''
    Filter mode tells you how you deal with oscillations in prediction
    If filter_mode=None, don't' do anything
    If filter_mode='First commit', then we insist that we commit to a change the first time it occurs
        For instnace, 01011 becomes 01111
        
    If filter_mode='Last commit', then we insist that we commit to a change the last time it occurs
        This is accomplished by using the function filter_oscillations with a thresh setting of 200, so
        any oscillation is filtered
    '''
    files, trajectories=visualize_PCA.Get_trajectory(data, PDB_files, traj)
    labels=Markov_prediction(trajectories, files, gmix, filter_mode=filter_mode)
    plt.figure()
    plt.plot(Get_times(files),labels)
    


def Get_times(files):
    times=[]
    for file in files:
        entries=file.split('.')
        times.append(float(entries[-1]))
    return times

def Plot_unfolding_trajectory(traj, labels, PDB_files, filter_osc=True, thresh=7, islands_only=True):
    f, l = visualize_PCA.Get_trajectory(labels, PDB_files, traj)
    if filter_osc==True: l=filter_oscillations(l, thresh=thresh, islands_only=islands_only)
    #if filter_osc==True: l=filter_oscillations(l, thresh=thresh)
    plt.figure()
    plt.plot(Get_times(f),l)
    
    unique=list(set(labels))
    if np.nan in unique: unique.remove(np.nan)
    plt.ylim((np.min(unique)-0.1, np.max(unique)+0.1))
    
def Cluster_by_residue_distance(thresh=8):
    [distances, PDB_files]=joblib.load('ADK_trunc130_tighter_temps/Residue_distances.dat')
    PDB_files, distances=visualize_PCA.sort_projections(PDB_files, distances)
    labels=[]
    for d in distances[:,1]:
        if d<thresh: labels.append(0)
        else: labels.append(1)
    pickle.dump([ labels, PDB_files], open("ADK_trunc130_tighter_temps/Sheet_folding_labels.dat", "wb"))
    return labels
    
    
    
def Combine_clusters(i,j,labels):
    """
    All clusters with label j are reassigned to cluster i, while all cluster labels above j are shifted down
    """
    if i>j:   #swap i and j
        i_old=i
        i=j
        j=i_old
    for n, l in enumerate(labels):
        if l==j: 
            labels[n]=i
        elif l>j:
            labels[n]=l-1
    return labels
    

def load_data(protein, include='psn', directories=['tighter_temps', 'unfolding'], n_pcs=10, pca_name='PCA_combined.dat', score_file='Substructure_scores.dat'):
    """
    You indicate what sort of data you want to include in this data matrix
    direcotry tells you where all this data is stored, should all be in same place
    
    
    d stands for data (energy, rmsd, and contacts) (NOT INCLUDED BY DEFAULT)
    p stands for PCA projections
    s stands for scores 
    n stands for native contacts
    o stands for non-native contacts
    
    n_pcs tells you how many principal components to include, if projetions is included
    
    Every type of data must have the same number of points..the data-types are concatenated into one big array whose vertical axis represents timepoints
    """

    directory=directories[0] #the main directory where everything is stored
    if 'p' in include:
        [ pca, projections, PDB_files]=visualize_PCA.load_PCA('{}_{}/{}'.format(protein, directory, pca_name))
        combined_data=projections[:,0:n_pcs]
    
    if 's' in include:
        scores_and_natives, PDB_files=load_scores(protein, directory=directory, score_file=score_file, elim_ID=False)
        scores=scores_and_natives[:, 0:-2] #omit last two columns which are natives and non-native contacts
        try:
            combined_data=np.hstack((combined_data, scores))
        except NameError: 
            combined_data=scores
            
    if 'n' in include:  #Generally, n_natives is the second to last column of the substructure scores matrix
        try: #in case we have already loaded the substructure scores
            natives=scores_and_natives[:,-2]
        except NameError: #in case we haven't
            natives, PDB_files=load_scores(protein, score_file=score_file)
            natives=natives[:,-2]
        natives=np.reshape(natives, (len(natives), 1))
        try: 
            combined_data=np.hstack((combined_data, natives))
        except NameError:
            combined_data=natives
    if 'o' in include:  #Generally, n_nonnatives is the last column of the substructure scores matrix
        try: #in case we have already loaded the substructure scores
            nonnatives=scores_and_natives[:,-1]
        except NameError: #in case we haven't
            nonnatives, PDB_files=load_scores(protein, score_file=score_file)
            nonnatives=nonnatives[:,-1]
        nonnatives=np.reshape(nonnatives, (len(nonnatives), 1))
        try: 
            combined_data=np.hstack((combined_data, nonnatives))
        except NameError:
            combined_data=nonnatives
    if 'd' in include:
        for n,direc in enumerate(directories):
            [x, temperatures, log_files]=visualize_PCA.load_Data('{}_{}/All_Data.dat'.format(protein, direc))
            if n==0:
                data=x
            else:
                data=np.vstack((data, x))
        try:
            combined_data=np.hstack((combined_data, data))
        except NameError:
            combined_data=data
                


    combined_data, PDB_files=Find_identical_trajectories.Eliminate_identical_trajectories(combined_data, PDB_files)
    return combined_data, PDB_files



def load_data_old(protein, include='dp', n_pcs=10, pca_name='PCA_combined.dat',dir_types=['replica_tighter_temps', 'unfolding'], score_file='Substructure_scores.dat'):
    """
    You indicate what sort of data you want to include in this data matrix
    d stands for data (energy, rmsd, and contacts)
    p stands for PCA projections
    s stands for scores
    
    n_pcs tells you how many principal components to include, if projetions is included
    
    It is assumed that the number of PC projections equals the total number of datapoints (contacts, energies, rmsd)
    """
    replica_tighter_temps='{}_tighter_temps'.format(protein)
    unfolding='{}_unfolding'.format(protein)
    replica=protein
    if 'd' in include:
        
        [data, temperatures, log_files]=visualize_PCA.load_Data('{}/All_Data.dat'.format(replica_tighter_temps))
        
        if 'unfolding' in dir_types:  
            [data2, unfolding_temps, unfolding_log_files]=visualize_PCA.load_Data('{}/All_Data.dat'.format(unfolding))
            
            if 'replica' in dir_types:
                [data3, rep_temps, rep_files]=visualize_PCA.load_Data('{}/All_Data.dat'.format(replica))
                data=np.concatenate((data, data2, data3), axis=0)
    
            else:
                data=np.concatenate((data, data2), axis=0) 
        if 'p' in include:
            [ pca, projections, PDB_files]=visualize_PCA.load_PCA('{}/{}'.format(replica_tighter_temps, pca_name))
        #[data3, rep_temps, rep_files]=visualize_PCA.load_Data('{}/All_Data.dat'.format(replica))
            combined_data=np.concatenate((data, projections[:,0:n_pcs]), axis=1)
        else:
            combined_data=data
            PDB_files=log_files+unfolding_log_files
    elif 'p' in include and 'd' not in include:
        [ pca, combined_data, PDB_files]=visualize_PCA.load_PCA('{}/PCA_combined.dat'.format(replica_tighter_temps))
        combined_data=combined_data[:,0:n_pcs]
    else: 
        print('Variable "include" must contain "d" and/or "p"!')
    
    if 's' in include:
        scores, PDB_files=load_scores(protein, score_file=score_file)
        combined_data=np.hstack((combined_data, scores))
    return combined_data, PDB_files





def load_scores(protein, directory='tighter_temps',score_file='Substructure_scores.dat', distinguish_traps=False, thresh=2, elim_ID=True, convert_to_binary=True, min_nonnative_substructures=2, gmix=False):
    """
    Thresh tells you how far residues in a substructure have to be relative to native distance, on average, for substructure to be considered formed
    if elim_ID is true, you eliminate identical trajectories
    """
    
    
    try: 
        Scores, N_nonnative_substructures, PDB_files, native_contacts, substructures= joblib.load('{}/{}/{}'.format(protein, directory, score_file))
        
        
        
        a, Scores=visualize_PCA.sort_projections(PDB_files, Scores)
        a, N_nonnative_substructures=visualize_PCA.sort_projections(PDB_files, N_nonnative_substructures)

        PDB_files=a
    except ValueError:
        Scores, PDB_files, native_contacts, substructures=joblib.load('{}/{}/{}'.format(protein, directory, score_file))
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
        


def Save(gmix, data, labels, PDB_files, protein, description='', directory='tighter_temps'):
    path='{}_{}/Gaussian_labels.dat'.format(protein, directory)
    pickle.dump([gmix, data, labels, description, PDB_files], open(path, "wb"))



def Combine_scores(path1, path2):
    """
    Loads unfolding score data in path1 and path2, and combines them into a single file by offsetting the trajectory nubmers in path2,
    then resaves as path1
    """
    Scores1, N_nonnative_substructures1, PDB_files1, native_contacts, substructures= joblib.load(path1)
    Scores2, N_nonnative_substructures2, PDB_files2, native_contacts, substructures= joblib.load(path2)
    
    
    traj_nums=[]  #x. file 'MARR/unfolding3/marr_1.000_399.99500000' will b given traj_num = 399 
    for f in PDB_files1:
        fields=f.split('_')
        traj_num=int(fields[-1].split('.')[0])
        traj_nums.append(traj_num)
        
    offset=np.max(traj_nums)+1  #to every traj number in PDB_files2, we add this value
    
    for i, F in enumerate(PDB_files2):
        fields=F.split('_')
        
        traj_num, time = fields[-1].split('.')
        traj_num=int(traj_num)+offset
        time=int(time)
        
        
        newstring=fields[0]
        for field in fields[1:-1]:
            newstring='{}_{}'.format(newstring, field)
        newstring='{}_{}.{}'.format(newstring, traj_num, time)
        PDB_files2[i]=newstring
    
    
    Scores=np.vstack((Scores1, Scores2))
    PDB_files=PDB_files1+PDB_files2
    N_nonnative_substructures= np.hstack((N_nonnative_substructures1, N_nonnative_substructures2))
    
    joblib.dump([Scores, N_nonnative_substructures, PDB_files, native_contacts, substructures], path1 )
    



#f=lookup_files('00001000', scores, PDB_files, traj='0.600_*')
#nonnatives_low=[n_nonnatives[i] for i in range(len(PDB_files)) if PDB_files[i] in f]
#
#f=lookup_files('00001000', scores, PDB_files, traj='0.900_*')
#nonnatives_high=[n_nonnatives[i] for i in range(len(PDB_files)) if PDB_files[i] in f]