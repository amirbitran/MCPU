# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 18:05:19 2017

@author: amirbitran
"""

import pickle
import Analyze_structures



#First, we analyze the full ADK protein, and we care about all the residues


native_file='ADK/replica/adk_0.100.0'
directory='ADK/replica_tight_spacing'
fileroot='adk'
min_temp=0.450
max_temp=0.915
n_temps=32
ignore=26
print(fileroot)

Full_ADK=Analyze_structures.Simulation_analysis(native_file, directory, fileroot, min_temp, max_temp, n_temps, ignore )
Full_ADK.get_native()
Full_ADK.analyze_trajectories()



pickle.dump(Full_ADK, open("{}/Full_ADK.dat".format(directory), "wb"), pickle.HIGHEST_PROTOCOL)  #HIGHEST PROTOCOL is needed to dump objects




#Now, we analyze ADK with the last 26 residues truncated

native_file='ADK_trunc26/replica/adk_t26_0.100.0'
directory='ADK_trunc26/replica_tight_spacing'
fileroot='adk_t26'
min_temp=0.450
max_temp=0.915
n_temps=32
ignore=0   #ignore here refers to how many residues we ignore. 26 residues have already been truncated from the protein, so we care about all the residues that are already there
print(fileroot)


ADK_trunc26=Analyze_structures.Simulation_analysis(native_file, directory, fileroot, min_temp, max_temp, n_temps, ignore )
ADK_trunc26.get_native()
ADK_trunc26.analyze_trajectories()


pickle.dump(ADK_trunc26, open("{}/ADK_trunc26.dat".format(directory), "wb"), pickle.HIGHEST_PROTOCOL)  #HIGHEST PROTOCOL is needed to dump objects



#Now, we analyze the full length DHFR



native_file='DHFR/replica/dhfr_0.100.0'
directory='DHFR/replica_tight_spacing'
fileroot='dhfr'
min_temp=0.850
max_temp=1.000
n_temps=31
ignore=10
print(fileroot)

Full_DHFR=Analyze_structures.Simulation_analysis(native_file, directory, fileroot, min_temp, max_temp, n_temps, ignore )
Full_DHFR.get_native()
Full_DHFR.analyze_trajectories()


pickle.dump(Full_DHFR, open("{}/Full_DHFR.dat".format(directory), "wb"), pickle.HIGHEST_PROTOCOL)  #HIGHEST PROTOCOL is needed to dump objects



#Finally, we analyze the DHFR protein with the last 10 residues truncated


native_file='DHFR_trunc10/replica/dhfr_t10_0.100.0'
directory='DHFR_trunc10/replica_tight_spacing'
fileroot='dhfr_t10'
min_temp=0.800
max_temp=0.955
n_temps=32
ignore=0
print(fileroot)

DHFR_trunc10=Analyze_structures.Simulation_analysis(native_file, directory, fileroot, min_temp, max_temp, n_temps, ignore )
DHFR_trunc10.get_native()
DHFR_trunc10.analyze_trajectories()




pickle.dump(DHFR_trunc10, open("{}/DHFR_trunc10.dat".format(directory), "wb"), pickle.HIGHEST_PROTOCOL)  #HIGHEST PROTOCOL is needed to dump objects
