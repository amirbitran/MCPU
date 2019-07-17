# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:08:42 2017

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
import Bio
from Bio import PDB
import pickle
import ClusterPCA

#directory='ADK_trunc26_tighter_temps/PCA_complete.dat'



def load_PCA(directory):
    """
    loads PCA results and sorts files/projections in order
    Also works for kinetic simulations
    """

    [pca, projections, mean, PDB_files]=joblib.load(directory)
    if type(PDB_files[0])==list:  #some of my PCA results had PDB_files saved as a list of lists
        PDB_files=[f for sublist in PDB_files for f in sublist]
    
    PDB_files, projections=sort_projections(PDB_files, projections)
    
    return pca, projections, PDB_files


def load_Data(path):
    """
    loads data and sorts files/data in order
    data is returned in much the same way as PCA results are in load_PCA...
    the first column is RMSDs, second is contacts, third is energies
    """
    [energies, contacts, rmsd, log_files, temperatures]=joblib.load(path)
    
    x, energies=sort_projections(log_files, energies)
    x, contacts=sort_projections(log_files, contacts)
    x, rmsd=sort_projections(log_files, rmsd)
    log_files, temperatures=sort_projections(log_files, temperatures)
    
    E=np.reshape(energies, (np.size(energies),1))
    C=np.reshape(contacts, (np.size(contacts),1))
    R=np.reshape(rmsd, (np.size(rmsd),1))
    
    if np.ndim(temperatures)==2:
        temperatures=[float(t[0]) for t in temperatures]
    
    data=np.concatenate((R,C,E),axis=1 )
    return data, temperatures, log_files






def load_Data_new(path, reverse_order=False):
    """
    loads data and sorts files/data in order
    Use this function for data that was read ever since I started doing umbrella sampling

    reverse_order=True means you want to sort in order of decreasing setpoints
    """
    [data, variables, log_files, temperatures, setpoints, times]=joblib.load(path)
   
    
    
    #Figure out how to re-sort data to be in order of increasing temperatures, and increasing setpoints within each temperature
    #ex. adk_T_0.400_10.log, adk_T_0.400_20.log, ... adk_T_0.400_200.log, adk_T_0.450_10.log, adk_T_0.450_20.log, etc...
    order=[]
    unique_temperatures=np.array(list(set(temperatures))) 
    unique_temperatures=unique_temperatures[np.argsort(unique_temperatures)] #sort the temperatures first
    for temp in unique_temperatures:
        indices_at_temp=np.array([t for t in range(len(temperatures)) if temperatures[t]==temp ])
        
        if reverse_order:
            order_at_this_temp=indices_at_temp[np.argsort(-1*np.array(setpoints)[indices_at_temp])]
        else:
            order_at_this_temp=indices_at_temp[np.argsort(np.array(setpoints)[indices_at_temp])]
        for o in order_at_this_temp:
            order.append(o)  
        
    #order=np.argsort(reporters)
    data=data[order, :, :]
    log_files=[log_files[x] for x in order]
    temperatures=np.array(temperatures)[order]
    setpoints=np.array(setpoints)[order]
    return data, temperatures, setpoints, log_files, times, variables


def Visualize_starting_files(path):
    data, temperatures, setpoints, log_files, times, variables=load_Data_new(path, reverse_order=True)
    natives=data[:,:,0]
    natives=natives.flatten()
    
    times_per_setpoint=int(len(natives)/len(setpoints))
    ramp=[s for s in setpoints for i in range(times_per_setpoint) ]
    
    plt.figure()
    plt.plot(natives, label='Native contacts')
    plt.plot(ramp, color='red', label='Setpoint')
    plt.title('Natives ramp', fontsize=30)
    plt.xlabel('time', fontsize=24)
    plt.ylabel('Number of native contacts', fontsize=24)
    plt.tick_params(labelsize=20)
    plt.legend(fontsize=24)


def Plot_substructure_unfolding(protein, traj, directory='unfolding', thresh=12):
    """
    Plot an unfolding trajectory by plotting the nubmer of formed substructures as a function of time
    Also you label the configuration whenever it changes
    """
    scores, PDB_files=ClusterPCA.load_scores(protein, directory=directory, thresh=thresh)
    files, X =Get_trajectory(scores, PDB_files, traj)
    times=ClusterPCA.Get_times(files)
    N=[]
    label_positions=[]
    labels=[]
    y_values=[]  #keep track to avoid collisions
    for i, x in enumerate(X):
        N.append(np.sum([int(j) for j in x]))
        if i==0:
            labels.append(x)
            y=N[i]+0.01
            label_positions.append((times[i], y))
            y_values.append(y)
            k=0
        else:
            if x!=X[i-1] and x not in labels: #label has changed
                labels.append(x)
                if N[i]>=N[i-1]:
                    tentative_y_value=N[i]+0.3
                    while tentative_y_value in y_values:
                        tentative_y_value=tentative_y_value+0.3
                    y_value=tentative_y_value
                else:
                    tentative_y_value=N[i]-0.3
                    if tentative_y_value in y_values:
                        y_value=tentative_y_value-0.3
                    else:
                        y_value=tentative_y_value
                y_values.append(y_value)
                label_positions.append((times[i], y_value))
                k+=1
     
    
    N=np.array(N)
    y_values=np.array(y_values)
    plt.figure()
    plt.plot(times, N)
    ymin=np.min((np.min(y_values), np.min(N)))-0.5
    ymax=np.max((np.max(y_values), np.max(N)))+0.5
    plt.ylim((ymin, ymax))
    
    plt.xlabel('time (MC steps)', fontsize=24)
    plt.ylabel('Number of formed substructures', fontsize=24)
    plt.tick_params(labelsize=20)
    for j in range(len(labels)):
        plt.annotate(labels[j], xy=label_positions[j], fontsize=20)
    
    
    


def Fix_ADK_Data(directory):
    """
    Since ADK is processed differently, you must apply this function to any All_Data file involving ADK before its' compatible with the rest of analysis
    """
    path='{}/All_Data.dat'.format(directory)    
    [energies, contacts, rmsd, temperatures]=joblib.load(path)
    temperature_strings=[]
    for T in temperatures:
        current_tempstring='{}'.format(T)
        while len(current_tempstring)<5: current_tempstring='{}0'.format(current_tempstring)
        temperature_strings.append(current_tempstring)
    
    if 'trunc' in path:
        fields=path.split('trunc')
        truncation=int(fields[1].split('_')[0])
        log_files=['adk_t{}_{}.log'.format(truncation,T) for T in temperature_strings]
    else:
        log_files=['adk_{}.log'.format(T) for T in temperature_strings]
    length=np.shape(energies)[1]
    freq=int(length/200) #down sample energies, contacts, rmsd
    
    energies=energies[:,0::freq]
    contacts=contacts[:,0::freq]
    rmsd=rmsd[:,0::freq]
    pickle.dump([energies, contacts, rmsd, log_files, temperatures],open("{}/All_Data.dat".format(directory), "wb"), -1)

def load_Data_raw(path):
    """
    Same as above, except that you now return energies, contacts, rmsd as separate variables
    """
    [energies, contacts, rmsd, log_files, temperatures]=joblib.load(path)
    x, energies=sort_projections(log_files, energies)
    x, contacts=sort_projections(log_files, contacts)
    x, rmsd=sort_projections(log_files, rmsd)
    log_files, temperatures=sort_projections(log_files, temperatures)
    
    
    if type(temperatures[0])!=float:
        temperatures=[float(t[0]) for t in temperatures]

    return energies, contacts, rmsd, log_files, temperatures



def Replica_2D_histogram(data, temperatures, eq=100):
    energies=data[:,eq:,1]
    contacts=data[:,eq:,0]
    
    
    colors=['r', 'b', 'g', 'y', 'k', 'm', 'c']
    #colormaps=['binary', 'pink' , 'spring', 'summer', 'autumn', 'winter', 'Wistia', 'cool', 'hot', 'copper', 'bone', 'afmhot' ]
    plt.figure()
    plt.xlim(np.min(contacts), np.max(contacts))
    plt.ylim(np.min(energies), np.max(energies))
    
    
    nodes_per_temp=int(len(temperatures)/len(np.unique(temperatures)))

    colors_used=[]
    for i in range(np.shape(energies)[0]):
        if i>0:
            Go=True
            while Go:
                color=np.random.choice(colors)
                if color!=colors_used[i-1]:
                    if i<=nodes_per_temp:
                        Go=False
                    else:
                        if color!=colors_used[i-nodes_per_temp]:
                            Go=False
            colors_used.append(color)
        else:
            color=np.random.choice(colors)
            colors_used.append(color)
            
        #cmap_index=int(np.mod(i, len(colormaps)))
        #color_index=int(np.mod(i, len(colors)))
        plt.scatter(contacts[i,:], energies[i,:], color=color)
        #plt.hist2d(contacts[i,:], energies[i,:], bins=50, cmap='binary')
        #plt.hist2d(contacts[i,:], energies[i,:], bins=50, cmap=colormaps[cmap_index])
        plt.xlabel('Contacts', fontsize=24)
        plt.ylabel('Energy', fontsize=24)
        plt.tick_params(labelsize=20)


def sort_projections(PDB_files, projections):
    sort=natsort.natsorted(zip(PDB_files, projections))
    PDB_files=[tuple[0] for tuple in sort]
    if type(projections)==list:
        projections=[tuple[1] for tuple in sort]
    else:
        projections=np.vstack([tuple[1] for tuple in sort])
    return PDB_files, projections
    
def ScatterPlot2D(projections, comp1=0, comp2=1, axis1='rmsd', axis2='contacts'):
    plt.figure()
    plt.scatter(projections[:,comp1], projections[:,comp2], label='All points')
    plt.xlabel('{}'.format(axis1), fontsize=25)
    plt.ylabel('{}'.format(axis2), fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=22, pad=2)
    plt.title('Data', fontsize=25)
    
def ScatterPlot3D(projections):
    fig=plt.figure()
    ax=fig.add_subPllot(111, projection='3d')
    ax.scatter(projections[:,0], projections[:,1], projections[:,2])
    ax.set_xlabel('Principal component 1')
    ax.set_ylabel('Principal component 2')
    ax.set_zlabel('Principal component 3')
    ax.set_title('Data projected onto principal components')

def num2str(num):
    string=str(num)
    while len(string)<5:
        string='{}0'.format(string)
    return string




def Get_trajectory(projections, PDB_files, traj_list):
    """
    Obtains PDB files that at some temperature or set of temperatures
    traj_list can either be a numpy array of temperature values (ex. np.arange),
    a single string corresponding to am temperature  (ex. '0.600'),
    a string with wildcards (ex. '0.6**'), or a list of strings, potentially with wildcards 
    (ex. ['0.6**', '0.7**'])
    
    Also works for kinetic simulations...just enter the trajectory number
    Ex. if you want files of the form adk_0.800_3._ _ _, just enter '3' for traj_list
    """
    if type(traj_list)==str: traj_list=[traj_list]
    if type(traj_list)==np.ndarray:
        traj_list=[num2str(t) for t in traj_list]
    #traj=[f for f in range(len(PDB_files)) for traj_num in traj_list if fnmatch.fnmatch( PDB_files[f], '*_{}*'.format(traj_num)) ]  #sequence of indices corresponding to current PDB file 
    traj=[f for f in range(len(PDB_files)) for traj_num in traj_list if fnmatch.fnmatch( PDB_files[f], '*{}*'.format(traj_num)) ]  #sequence of indices corresponding to current PDB file 
    traj_files=[PDB_files[f] for f in traj]
    traj_coords=np.array([projections[f] for f in traj])
    return traj_files, traj_coords
    
    
#plt.figure()
  
def ScatterTrajectory3D(projections, PDB_files, traj_num):
    traj_files,f, traj_coords = Get_trajectory(projections,PDB_files, traj_num)
    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    ax.scatter(traj_coords[:,0], traj_coords[:,1], traj_coords[:,2])
    #plt.xlim(-400,1000)
    #plt.ylim(-400,500)
    ax.set_xlabel('Principal component 1')
    ax.set_ylabel('Principal component 2')
    ax.set_zlabel('Principal component 3')
    #plt.title('Scatter plot with trajectory {}'.format(traj_num))
    for n, file in enumerate(traj_files): print('{} \n {} \n'.format(file, (traj_coords[n,0], traj_coords[n,1])))
    #return traj_files, traj_coords
  
def ScatterTrajectory(projections, PDB_files, traj_num, comp1=0, comp2=1):
    plt.figure()
    traj_files, traj_coords = Get_trajectory(projections,PDB_files, traj_num)
    plt.scatter(traj_coords[:,comp1], traj_coords[:,comp2])
    #plt.xlim(-400,1000)
    #plt.ylim(-400,500)
    plt.xlabel('Principal component 1')
    plt.ylabel('Principal component 2')
    #plt.title('Scatter plot with trajectory {}'.format(traj_num))
    for n, file in enumerate(traj_files): print('{} \n {} \n'.format(file, (traj_coords[n,0], traj_coords[n,1])))
        
def PlotTrajectory(projections, PDB_files,traj_num, background='On', dim1=0, dim2=1):
    
    if background=='On': 
        ScatterPlot2D(projections)
    else:
        plt.figure()
        
    traj_files, traj_coords = Get_trajectory(projections,PDB_files, traj_num)
    plt.plot(traj_coords[:,dim1], traj_coords[:,dim2])
    N=np.shape(traj_coords)[0]
    for i in range(N-1):
        plt.plot(traj_coords[i:i+2,dim1], traj_coords[i:i+2,dim2], color=plt.cm.jet(i/N))
    #plt.xlim(-400,1000)
    #plt.ylim(-400,500)
    plt.xlabel('Principal component {}'.format(dim1))
    plt.ylabel('Principal component {}'.format(dim2))
    plt.title('Scatter plot with trajectory {}'.format(traj_num))
    for n, file in enumerate(traj_files): print('{} \n {} \n'.format(file, (traj_coords[n,dim1], traj_coords[n,dim2])))
    #return traj_files, traj_coords

def PlotTrajectory3D(projections, PDB_files,traj_num, dim1=0, dim2=1, dim3=2):
    fig=plt.figure()
    traj_files, traj_coords = Get_trajectory(projections,PDB_files, traj_num)
    ax=fig.add_subplot(111, projection='3d')
    ax.plot(traj_coords[:,dim1], traj_coords[:,dim2])
    N=np.shape(traj_coords)[0]
    for i in range(N-1):
        ax.plot(traj_coords[i:i+2,dim1], traj_coords[i:i+2,dim2], traj_coords[i:i+2,dim3], color=plt.cm.jet(i/N))
    #plt.xlim(-400,1000)
    #plt.ylim(-400,500)
    plt.xlabel('Principal component {}'.format(dim1))
    plt.ylabel('Principal component {}'.format(dim2))
    plt.title('Scatter plot with trajectory {}'.format(traj_num))
    for n, file in enumerate(traj_files): print('{} \n {} \n'.format(file, (traj_coords[n,dim1], traj_coords[n,dim2])))
    #return traj_files, traj_coords



def PlotTrajectoryTimecourse(projections, PDB_files, traj_num,  component_to_plot):
    traj_files, traj_coords =   Get_trajectory(projections, PDB_files, traj_num)  
    plt.figure()
    plt.plot(traj_coords[:,component_to_plot])
    plt.xlabel('MC step')
    plt.title('Position along principal component {}'.format(component_to_plot+1))



def Visualize_PCs_magnitude(pca):
    """
    plots the magnitude of the PCA component at each residue
    """
    n_residues=int(np.shape(pca.components_)[1]/3)
    
    magnitudes=np.zeros((np.shape(pca.components_)[0], n_residues))
    plt.figure()
    
    for j in range(np.shape(pca.components_)[0]):
        for res_n, i in enumerate(np.arange(0, np.shape(pca.components_)[1], 3 )):
            magnitudes[j, res_n]= np.sqrt(np.dot(pca.components_[j, i:i+3], pca.components_[j, i:i+3]))
        plt.plot(np.arange(1,n_residues+1,1), magnitudes[j,:], label='Component {}'.format(j+1))
    plt.xlabel('Residue number')
    plt.xlim(0,n_residues)
    plt.legend()
    plt.title('Magnitude of principal components at residues')
            
def Visualize_PCs(pca, reference_file):
    """
    plots the radial component of the three PCA components at each residue
    to get a sense of both magnitude and direction, we normalize the radial vector but do not normalize the three PCA components at each residue
    Radial is defined relative to the mean coordinate of the reference structure, assumed to be the native structure
    Positive values means the component points radially outwards, while negative means inwards

    """
    
    #First, load the reference strcture and obtain the coordinates for the center of the protein
    reference=Bio.PDB.PDBParser().get_structure('Protein',reference_file).get_residues()
    ref_structure=[]
    for res in reference: ref_structure.append(res['CA'])
    ref_coords=np.array([list(atom.get_coord()) for atom in ref_structure])


    protein_center=np.mean(ref_coords, axis=0)    
    
    #now we 
    n_residues=np.shape(ref_coords)[0]
    
    radial_components=np.zeros((np.shape(pca.components_)[0], n_residues))
    plt.figure()
    
    for j in range(np.shape(pca.components_)[0]):   #loop throuhg the different PCA components
        for res_n, i in enumerate(np.arange(0, np.shape(pca.components_)[1], 3 )):  #lloop through the residues
            radial_vector=ref_coords[res_n,:]-protein_center
            radial_vector=radial_vector/np.sqrt(np.dot(radial_vector, radial_vector))  #normalize
            
            pca_vector=pca.components_[j, i:i+3]
            radial_components[j, res_n]=np.dot(radial_vector, pca_vector)
        plt.plot(np.arange(1,n_residues+1,1), radial_components[j,:], label='Component {}'.format(j+1))
    plt.xlabel('Residue number')
    plt.xlim(0,n_residues)
    plt.legend()
    plt.title('Radial projection of principal components at residues')
            


def Create_PC_Decoy(ref_directory, ref_structure, component):
    """
    Prepares a PDB structure that can be loaded in PYMOL along with a reference strucutre to visualize  principal components
    The file loads a reference structure (ex. a native strucutre) found in ref_directory/ref_structure, as well as the PCA results stored in ref_directory
    It then transforms the atomic coordinates of ref_structure by adding, to each residue, the value of some principal component at that residue
    The principal component number is specified by component (ex. component=0)
    
    Note: ref_structure should not have .pdb at the end
    """
    reference=Bio.PDB.PDBParser().get_structure('Protein','{}/{}.pdb'.format(ref_directory, ref_structure))
    transformed=reference   #A structure that we will transform by adding the eigenvector 
    [pca, projections, PDB_files]=load_PCA('{}/PCA_complete.dat'.format(ref_directory))
    for model in transformed:
        for chain in model:
                for n, residue in enumerate(chain):
                    for atom in residue:
                        pca_vector=pca.components_[component, 3*n:3*n+3]   #PCA vector at this residue, has x y and z coordinates
                        atom.set_coord(atom.get_coord()+50*pca_vector)   #We transform each atom by adding the PCA vector at that residue (up to a factor of 50 to make the change more notable)
                        #atom.set_coord([0,0,0])  #The transformation is simply setting all atomic coordinates equal to (0,0,0)
    
    io=PDB.PDBIO()
    io.set_structure(transformed)    
    io.save('{}/{}_comp{}.pdb'.format(ref_directory,ref_structure, component )) 



########################################################



    




#reference_file='ADK/1ake.pdb'

#[ pca, projections, PDB_files]=load_PCA(directory)
#ScatterPlot2D(projections)

#[pca2, projections2, mean2, PDB_files2]=joblib.load('ADK_unfolding/Secondary_PCA.dat')  #we plot the secondary PCA results done on 

