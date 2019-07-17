import numpy as np
import os
import joblib
import itertools
import joblib
import Analyze_structures
import LoopCluster
import glob
import copy as cp
import natsort
import fnmatch
import copy as cp
import ClusterSubstructures
import ClusterPCA


protein='MARR_trunc44'#ex. MARR or MARR_trunc44
fileroot='marr_t44' #ex. marr or marr_t44

min_step = 0 #only read steps that are greater than this number (inclusive)
max_step=np.inf  #Only read steps that are less than this number
d_cutoff=7.8   #6.5 for FABG, DHFR, 1igd, and CMK, 7.8 for MARR , 8 for HEMK, 5.8 for GCVH, 7 for RSME, 5.5 for coarse RSME
min_seq_separation = 8
directory='MultiUmbrella7'#what directories we analyze
thresh=2.1  #distance thresh for score assignment


starting_config = '010000'  #look at snapshots assigned to this config
sub_to_form = 2 #Count nonnatives that need to be broken in order to form this substructure


#Substructure_file = '1igd/MultiUmbrella7/Substructure_scores.dat'
#Substructure_file='FABG/MultiUmbrella5/Substructure_scores.dat'
#Substructure_file='ADK/MultiUmbrella4/Substructure_scores.dat'
#Substructure_file='CMK/unfolding3/Substructure_scores.dat'
#Substructure_file='DHFR/MultiUmbrella/Substructure_scores.dat'
#Substructure_file = 'HEMK/files/Substructure_scores.dat'
#Substructure_file = 'GCVH/MultiUmbrella/Substructure_scores.dat'
Substructure_file = 'MARR/MultiUmbrella3/Substructure_scores.dat'
#Substructure_file = 'MARR_hairpin/MultiUmbrella2/Substructure_scores.dat'
#Substructure_file = 'RSME/MultiUmbrella3/Substructure_scores.dat'  #MultiUmbrella3 has the more coarse grained parameters, MultiUmbrella has the more detailed ones
#Substructure_file = None



alphabet = 'abcdefghijklmnopqrstuvwxyz'
#string = str([alphabet[i] for i, bit in enumerate(starting_config) if bit ==1])
output_filename='Nonnatives_impeding_{}.dat'.format(alphabet[sub_to_form])



native_file='{}/files/{}_0.100_Emin.pdb'.format(protein, fileroot)
temperatures=[['*.***']]




scores,PDB_files=ClusterPCA.load_scores(protein, distinguish_traps=False, min_nonnative_substructures=2, directory=directory, thresh=thresh)
unused1, unused2, unused3, unusd4, native_substructures = joblib.load(Substructure_file)


coords, resis=Analyze_structures.read_PDB(native_file, 'CA')
native_contacts=Analyze_structures.compute_contacts_matrix(coords, mode='binary', thresh = d_cutoff, min_seq_separation=min_seq_separation)

length = np.shape(native_contacts)[0]
if length < np.shape(native_substructures)[0]:
	native_substructures = native_substructures[0:length, 0:length, :]  #only keep whatever portions of substructures are acutally found in the current protein
	native_contacts=native_contacts[0:length, 0:length] #same for native contacts
substructure = native_substructures[:,:,sub_to_form]


def Count_nonnatives_to_break(snapshot, native_contacts, substructure, d_cutoff, min_seq_separation):
	"""
	Given a snapshot, finds non-native contacts involving residues that participate
	in substructure of interest
	Substructure should be an NxN array (N is number of amino acids) with 1's at contacts
	involved in substructure
	"""
	coords, resis=Analyze_structures.read_PDB(snapshot, 'CA')
	contacts=Analyze_structures.compute_contacts_matrix(coords, mode='binary', thresh = d_cutoff, min_seq_separation=min_seq_separation)
	inverse_native_contacts = -(native_contacts - np.ones(np.shape(native_contacts)))
	nonnatives = np.multiply(contacts, inverse_native_contacts) #nonnative contacts in snapshot
	indices_to_filter = np.where(substructure)
	Filter = np.zeros(np.shape(substructure))
	Filter[indices_to_filter[0],:] = 1
	Filter[indices_to_filter[1],:] = 1
	Filter[np.where(Filter)]=1  
	Filter = np.maximum(Filter, Filter.transpose())
	return np.sum(np.multiply(Filter, nonnatives))




N_nonnatives = np.nan*np.ones(len(PDB_files))
indices = np.array([i for i in range(len(PDB_files)) if scores[i]==starting_config])

print('{} files total'.format(len(indices)))
for n, ind in enumerate(indices):
	if np.mod(n, 500) ==0: print('{} files completed'.format(n))
	snapshot = PDB_files[ind]
	N_nonnatives[ind] = Count_nonnatives_to_break(snapshot, native_contacts, substructure, d_cutoff, min_seq_separation)


joblib.dump([PDB_files, N_nonnatives, native_substructures, starting_config, sub_to_form], open('{}/{}/{}'.format(protein, directory, output_filename), "wb"), compress=3)




	