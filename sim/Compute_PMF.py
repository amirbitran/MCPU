# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 12:29:37 2017

@author: amirbitran

Computes potential of mean force (free energy) for protein as a function of fraction of native contacts formed at a chosen temperature,
using energies from replica exchange simulations at many temperatures

Usesn the script parallel-tempering-2dpmf.py, stored in /anaconda/lib/python3.5/site-packages/pymbar/testsystems

REFERENCES
[1] Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states.
J. Chem. Phys. 129:124105, 2008
http://dx.doi.org/10.1063/1.2978177
"""


import numpy as np
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib
import matplotlib.pyplot as plt



##############################################################################
# First, read the file
##############################################################################

directory='DHFR/replica_tight_spacing'
protein='dhfr'
min_filenum=0.850
max_filenum=1.000
n_files=31
n_truncated=0
n=10000  #how many samples from contacts and energies do we want? The fewer, the faster this thing runs


def read_file(directory,min_file, max_file, number_of_files, protein, truncated):
    
    filerange=np.linspace(min_file, max_file, number_of_files)
    energies=[]  #rows will correspond to replicas, columns to samples within a replica
    contacts=[]
    temperatures=[float(str(x)) for x in filerange]  #Really stupid step I have to do, because np.linspace  does silly things like output 0.15 as 0.1499999. FOr some reason, converting to a string and back fixes this   
    lens=[]  #length of the datasets: will be the same for all temperatures if all simulations are complete, but may vary if simulations are truncated
    for k, filenumber in enumerate(filerange):
        #filenumber=np.around(filenumber, decimals=1)
        print("Reading file {}".format(filenumber))        
        energies.append([])
        contacts.append([])
        filenumber=str(filenumber)
        while len(filenumber)<5:  #add zeros to the end of the filename until it has 3 zeros
            filenumber='{}0'.format(filenumber)
        n=0

        if n_truncated>0:
            openfile=open('{}/{}_t{}_{}.log'.format(directory, protein, n_truncated,filenumber))
        else:
            openfile=open('{}/{}_{}.log'.format(directory, protein, filenumber))
        #rmsd=[]
        for line in openfile.readlines():
            line=line.rstrip('\n')
            if len(line)>0:
                entries=line.split()
                if entries[0]=='STEP':
                    energies[k].append(float(entries[2]))
                    contacts[k].append(float(entries[3]))
            
                    #rmsd.append(float(entries[4]))
                    #print(k)                    
                    n+=1
        lens.append(len(contacts[k]))
        
        #contacts[k]=list(np.array(contacts[k])/max(contacts[k])) #temporarily convert to NP array to make it faster to normalize
        
        #contacts[k]=[i/max(contacts[k]) for i in contacts[k]]  #normalize so that gives  a fraction of native contacts
        
    contacts=np.array([x[0:min(lens)] for x in contacts])  #convert these lists to arrays, but in case the datasets are not all the same lenght, use the shortest length
    #contacts=contacts/np.max(contacts, axis=1)[:,None]      #normalize contactsto give a fraction of native contacts
    energies=np.array([x[0:min(lens)] for x in energies])
    return energies,contacts, temperatures  
                    
if os.path.exists("{}/All_Data.dat".format(directory)):
    print("Opening data")
    energies, contacts, temperatures=joblib.load("{}/All_Data.dat".format(directory))
else:

    print('Reading files')        
    energies, contacts, temperatures=read_file(directory, min_filenum,max_filenum,n_files, protein, n_truncated)

    print("Saving data")
    pickle.dump([energies, contacts, temperatures],open("{}/All_Data.dat".format(directory), "wb"), -1)
    


#minlen=min([len(x) for x in energies])
#X=[x[0:minlen] for x in contacts]
##############################################################################
# Now, downsample to obtain independent samples
##############################################################################

#contacts_reduced=contacts[:,-2701:-1]  #contacts is a huge matrix, so for testing purposes, we will work with a truncated version

def downsample(contacts,energies,n):
    """ 
    For now, this artificially downsamples the list of contacts and energies, by returning n evenly spaced samples
    This is just to speed up the process
    ToDo: Use legitimate statistical downsampling
    """


    indices=[int(x) for x in np.linspace(len(contacts[0])/2,len(contacts[0]),n) if x<len(contacts[0])]
    #indices=[int(x) for x in np.linspace(len(contacts[0]),len(contacts[0]),n) if x<len(contacts[0])]
#starting halfway through contacts matrix so that we ignore that equilibration time
    contacts_reduced=np.array([np.array(N)[indices] for N in contacts])
    energies_reduced=np.array([np.array(N)[indices] for N in energies])
    return contacts_reduced, energies_reduced  #These are numpy arrays

    #for n,x in enumerate(contacts_reduced):  #normalize to 
        #contacts_reduced[n,:]=x/np.max(x)
contacts_reduced, energies_reduced = downsample(contacts, energies, n)
   

print("Computing statistical inefficiencies...")
#g=timeseries.statisticalInefficiencyMultiple(contacts_reduced[:,5000:9999])

#snippet to play around with statistical inefficiency
#def experiment(o):
    #o=11
#    g=timeseries.statisticalInefficiencyMultiple(contacts_reduced[o,:])
#    plt.plot(contacts_reduced[o,:])
#    print(g)
#Continue for now
#ToDo: FIgure out how to obtain statistical inefficiency values that are much shorter than length of simulation
#This seemst o work much better if I omit simulations that have a relaxation time that is comperable to simulation time (ex. T=0.1) 

#ToDo: Downsample contacts and energies accordingly


#===================================================================================================
# Generate a list of indices of all configurations in kn-indexing
#===================================================================================================

#ToDo: Restructure by defining some of these constants earlier
# Create a list of indices of all configurations in kn-indexing.
N_k=[len(x) for x in contacts_reduced ]
N_max=max(N_k)  #note: this might not be correct once you downsample
K=len(temperatures)

def mask(N_k, N_max,K):
    mask_kn = np.zeros([K,N_max], dtype=np.bool)
    for k in range(0,K):
        mask_kn[k,0:N_k[k]] = True
# Create a list from this mask.
    indices = np.where(mask_kn)
    return mask_kn, indices
mask_kn, indices=mask(N_k, N_max, K)


#===================================================================================================
# Compute reduced potential energy of all snapshots at all temperatures
#===================================================================================================

print("Computing reduced potential energies...")

beta_k=[1/t for t in temperatures]  #Checked from simulations that temperatures are values of kBT
#Ex. a temperature of 0.1 means kBT=0.1

def reduced_potential_energies(energies_reduced, beta_k,K, N_max):
    u_kln = np.zeros([K,K,N_max], np.float32) # u_kln[k,l,n] is reduced potential energy of trajectory segment n of temperature k evaluated at temperature l
#pages correspond to different time points, rows correspond to different temperatures at whcih the time-series were taken, columsn correspond to temperatures at which energies shoudl be evaluated. Diagonal on any given page corresponds to original energies evaluated at their original temperatures at that time point 
    for k in range(K):
        for l in range(K):
            u_kln[k,l,0:N_k[k]] = beta_k[l] * energies_reduced[k,0:N_k[k]]
    return u_kln
u_kln=reduced_potential_energies(energies_reduced, beta_k,K, N_max)
      
    

#===================================================================================================
# Bin contact values into histogram bins for PMF calculation
#===================================================================================================

# Here, we bin the contacts samples into bins in a 2D histogram.
# We assign indices 0...(nbins-1) to the bins
# All bins must have at least one sample in them.
# This strategy scales to an arbitrary number of dimensions.

print("Binning energies...")
# Determine contacts bin size 
energies_min=-755
#contacts_max=1400
energies_max=-64
nbins_total=140


def bin(X, xmin, xmax, nbins_total):
    """ 
    Splits the observable X between values xmin and xmax into nbins_total bins    
    """
    
    
    dx = (xmax - xmin) / float(nbins_total)
    # Assign contacts bins
    bin_kn = np.zeros([K,N_max], np.int16) # bin_kn[k,n] is the index of which histogram bin does sample n from temperature index k belong to
    nbins = 0
    bin_counts = list()
    bin_centers = list() # bin_centers[i] gives the center of bin i
    for i in range(nbins_total):
    
      # Determine contacts value of bin center.
      xc = xmin + dx * (i + 0.5)
     
    #HERE!!!
      # Determine which configurations lie in this bin.
      in_bin = (xc-dx/2 <= X[indices]) & (X[indices] < xc+dx/2) 
    
      # Count number of configurations in this bin.
      bin_count = in_bin.sum()
    
      # Generate list of indices in bin.
      indices_in_bin = (indices[0][in_bin], indices[1][in_bin])
    
      if (bin_count > 0):
         # store bin (phi,psi)
         bin_centers.append( xc )
         bin_counts.append( bin_count )
    
         # assign these conformations to the bin index
         bin_kn[indices_in_bin] = nbins
    
         # increment number of bins
         nbins += 1
    
    print("%d bins were populated:" % nbins)
    for i in range(nbins):
       print("bin {} ({}) {} conformations".format(i, bin_centers[i], bin_counts[i]))
    return bin_kn, nbins, bin_centers

bin_kn, nbins, bin_centers=bin(energies, energies_min, energies_max,nbins_total)



binning_method="Binned based onenergies"
#print("Saving reduced data")



#===================================================================================================
# Compute PMF at the desired temperature.
#===================================================================================================



def test(min_temp, max_temp, temperatures, u_kln, N_k):
    temperatures=np.round(np.array(temperatures), decimals=2)
    temps_of_interest=np.linspace(min_temp, max_temp, round((max_temp-min_temp)/0.05)+1)
    temps_of_interest=np.round(temps_of_interest, decimals=2)
    indices=np.where(np.in1d(temperatures, temps_of_interest))[0]
    mbar=pymbar.MBAR(u_kln[indices[0]:indices[-1]+1, indices[0]:indices[-1]+1,:],N_k[indices[0]:indices[-1]+1], verbose=True )
    return mbar
    
     #mbar = pymbar.MBAR(u_kln[(len(u_kln)-n_temps_to_keep):len(u_kln),(len(u_kln)-n_temps_to_keep):len(u_kln),:], N_k[(len(N_k)-n_temps_to_keep):len(N_k)], verbose=True)


#t1=0.1
#t2=0.15

#temps_of_interest=np.array([t1, t2])
#temperatures=np.round(np.array(temperatures), decimals=2)
#indices=np.where(np.in1d(temperatures, temps_of_interest))[0]

#exp_diffs=np.exp(u_kln[indices[0], indices[0],:]-u_kln[indices[0], indices[1],:])

#(1/t1-1/t2)*energies_reduced[indices[0],:]


def compute_PMF(energies_reduced, contacts_reduced, temperatures, bin_kn, nbins, u_kln, N_k):
    """ Loops through every temperature and computes PMF
    u_kln are energies normalized by kT    
    """    
    print("Computing potential of mean force...")
    
    
    #energies_reduced, contacts_reduced, temperatures, binning_method, bin_kn, nbins, u_kln, N_k=joblib.load("{}/Reduced_Data.dat".format(directory))
    
    
    # Initialize MBAR
    mbar = pymbar.MBAR(u_kln, N_k, verbose=True)
    free_energies={}
    
    
    # Compute reduced potential energies at all temperatures
    for target_temperature in temperatures:
    #target_temperature=1.0
        target_temperature=float(target_temperature)        
        print("Computing PMF for temperature {}".format(target_temperature))        
        target_beta = 1.0 / ( target_temperature)
        u_kn = target_beta * energies_reduced
        # Compute PMF at this temperature, returning dimensionless free energies and uncertainties.
        # f_i[i] is the dimensionless free energy of bin i (in kT) at the temperature of interest
        # df_i[i,j] is an estimate of the covariance in the estimate of (f_i[i] - f_j[j], with reference
        # the lowest free energy state.
        (f_i, df_i) = mbar.computePMF(u_kn, bin_kn, nbins, uncertainties='from-normalization')
        f_i=f_i+np.log(np.sum(np.exp(-f_i)))   #ensures PMF's are normalized
        free_energies["{}".format(target_temperature)]=f_i
        
        
    
    # Show free energy and uncertainty of each occupied bin relative to lowest free energy
    #print("2D PMF")
    #print("")
    #print("%8s %6s %8s %10s %10s" % ('bin', 'contacts_reduced', 'N', 'f', 'df'))
    
    #for i in range(nbins):
    #   print('%8d %6.1f %8d %10.3f %10.3f' % (i, bin_centers[i], bin_counts[i], f_i[i], df_i[i]))
    #plt.plot(f_i)
    return free_energies


free_energies=compute_PMF(energies_reduced, contacts_reduced, temperatures, bin_kn, nbins, u_kln, N_k)
pickle.dump([energies_reduced, contacts_reduced, temperatures, binning_method, bin_kn, nbins, bin_centers, u_kln, N_k, free_energies],open("{}/Free_energies.dat".format(directory), "wb"), -1)

import matplotlib.pyplot as plt
import joblib

def plot_free_energies(temps_to_plot, directory, protein=None):
    if protein==None: protein=directory
    if temps_to_plot is not list: temps_to_plot=[temps_to_plot]
    energies_reduced, contacts_reduced, temperatures, binning_method, bin_kn, nbins, bin_centers, u_kln, N_k, free_energies=joblib.load("{}/Free_energies.dat".format(directory))
    #plt.figure(1)
    for temp in temps_to_plot:
        temp=temperatures[temperatures.index(temp)] #want to make sure temp has the same number of decimal points as in the list temperatures, since this is required to access it from dict
        plt.plot(bin_centers, free_energies[str(temp)], label=str(temp) )
    #plt.xlim(0,1)
    plt.ylim(0,25)
    plt.legend()
    plt.title("Potential of mean force, {}".format(protein))
    plt.xlabel("Fraction of native contacts formed")
        
        
    

