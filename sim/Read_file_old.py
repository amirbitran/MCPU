
import numpy as np
#import pymbar # for MBAR analysis
#from pymbar import timeseries # for timeseries analysis
import os
import pickle
#import os.path
import joblib



##############################################################################
# First, read the file
##############################################################################

directory='CRP_trunc110/unfolding'
protein='crp_t110'  #if protein is truncated, this should be reflected...ex. adk_t26
#min_temp=0.570
#max_temp=1.095
#spacing=0.015
temprange=[0.900, 1.000, 1.250, 1.750, 1.100, 1.200]
files_per_temp=9

file_numbers=np.arange(0,36,1) #This varible only matters if you have more than 1 file per temp
#In that case, your filenames will be of the form adk_0.800_1.pdb
#This array tells you what are the numbers that should come after the temperatures
#This array's length needs to equal (max_temp - min_temp)*files_per_temp/spacing
#that is, every file needs to have a number
#Again, this variable only plays any role if files_per_temp is not 1

def read_file(directory,min_temp, max_temp, spacing, protein, files_per_temp, file_numbers):
    temprange=np.arange(min_temp, max_temp+0.00001, spacing)   #+0.00001 is used to ensure that upper bound is included in range
    energies=[]  #rows will correspond to replicas, columns to samples within a replica
    contacts=[]

    temprange=[float(str(x)) for x in temprange]  #Really stupid step I have to do, because np.linspace  does silly things like output 0.15 as 0.1499999. FOr some reason, converting to a string and back fixes this   
    lens=[]  #length of the datasets: will be the same for all temperatures if all simulations are complete, but may vary if simulations are truncated
    rmsd=[]
    
    temperatures=[]  # temperature of each file, we fil this in as we go
    
    filecounter=0  #keeps track of which file you are currenatly reading
    
    for k, temp in enumerate(temprange):
        temp=str(temp)
        while len(temp)<5:  #add zeros to the end of the filename until it has 3 zeros
            temp='{}0'.format(temp)
        for j in range(files_per_temp):
            if files_per_temp>1:
                filenumber=file_numbers[filecounter]
                suffix='T_{}_{}.log'.format(temp, filenumber)
                temperatures.append('{}_{}'.format(temp, filenumber))
            else:
                suffix='{}.log'.format(temp)
                temperatures.append(float(temp))
            filename='{}/{}_{}'.format(directory, protein, suffix)
            print("Reading file {}".format(filename))      
            openfile=open(filename)
            
            energies.append([])
            contacts.append([])
            rmsd.append([])
            #temperatures.append(float(temp))
            for line in openfile.readlines():
                line=line.rstrip('\n')
                if len(line)>0:
                    entries=line.split()
                    if entries[0]=='STEP':
                        energies[filecounter].append(float(entries[2]))
                        contacts[filecounter].append(float(entries[3]))
                        rmsd[filecounter].append(float(entries[4]))
            lens.append(len(contacts[filecounter]))  
            filecounter+=1     

    contacts=np.array([x[0:min(lens)] for x in contacts])  #convert these lists to arrays, but in case the datasets are not all the same lenght, use the shortest length
    energies=np.array([x[0:min(lens)] for x in energies])
    rmsd=np.array([x[0:min(lens)] for x in rmsd])
    return energies,contacts, rmsd, temperatures  
                    
print('Reading files')        
energies, contacts, rmsd,temperatures=read_file(directory, min_temp ,max_temp, spacing, protein, files_per_temp, file_numbers )

print("Saving data")
pickle.dump([energies, contacts, rmsd, temperatures],open("{}/All_Data.dat".format(directory), "wb"), -1)
