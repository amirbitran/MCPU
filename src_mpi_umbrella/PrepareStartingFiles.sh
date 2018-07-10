#!/bin/bash
#First input is contacts max, second is contacts spacing, third is number of values to sweep over
#Fourth is path for input PDB file (Note that '/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/' is implicit)
#Fifth is path to output PDB files that are to be prepared (root above is implicit as well)...note that this should include PDB fileroot, like adk
#Sixth is the root of the protein, ex. 'adk'

#An example instance of running this would be . ./PrepareStartingFiles.sh 300 10 30 'ADK/files/adk_0.100_Emin.pdb' 'ADK/starting_files' 'adk'


((Contacts_max=$1))
((Contacts_spacing=$2))
((N=$3))
((Contacts_min=Contacts_max-N*Contacts_spacing))
protein_root=$6
root='/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/'
input_path="$root$4"
output_path="$root$5"

temperature="0.650"


#We use the awk command to edit files...the OFS tells it to use the three tabs as the delimiter
awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="NATIVE_FILE"){$2=var}} {print}' cfg_prepare>temp && mv temp cfg_prepare  #Edit native file
awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' cfg_prepare>temp && mv temp cfg_prepare  #Edit structure file
awk -v var="$output_path/$protein_root" '{OFS="\t\t\t"}{if ($1=="PDB_OUT_FILE"){$2=var}} {print}' cfg_prepare>temp && mv temp cfg_prepare  #Edit structure file
awk '{OFS="\t\t\t"}{if ($1=="NATIVE_DIRECTORY"){$2="None"}} {print}' cfg_prepare>temp && mv temp cfg_prepare
awk '{OFS="\t\t\t"}{if ($1=="MC_STEPS"){$2=100001}} {print}' cfg_prepare>temp && mv temp cfg_prepare
awk '{OFS="\t\t\t"}{if ($1=="MC_REPLICA_STEPS"){$2=100000000000000000000}} {print}' cfg_prepare>temp && mv temp cfg_prepare #No replica 
awk '{OFS="\t\t\t"}{if ($1=="MC_PDB_PRINT_STEPS"){$2=100000}} {print}' cfg_prepare>temp && mv temp cfg_prepare
awk '{OFS="\t\t\t"}{if ($1=="K_BIAS"){$2=0.02}} {print}' cfg_prepare>temp && mv temp cfg_prepare
awk '{OFS="\t\t\t"}{if ($1=="NODES_PER_TEMP"){$2=1}} {print}' cfg_prepare>temp && mv temp cfg_prepare
awk -v var="$temperature" '{OFS="\t\t\t"}{if ($1=="MC_TEMP_MIN"){$2=var}} {print}' cfg_prepare>temp && mv temp cfg_prepare  #Edit structure file





for ((i=0; i<N; i++)); do
	#Compute current centerpoint value
	((Contacts=Contacts_max-Contacts_spacing*i))
	
	#Edit cfg file so that simulation uses this centerpoint
	
	awk -v var="$Contacts" '{OFS="\t"}{if ($1=="NUMBER_OF_CONTACTS_MAX"){$2=var}} {print}' cfg_prepare>temp && mv temp cfg_prepare  #Edit native file
	
	#Run the simulation
	echo "Running simulation at Contacts=$Contacts"
	mpiexec -n 1 fold_potential_mpi cfg_prepare
	
	#now, we edit the cfg file so that the next simulation uses the output of the current one to start
	last_file_from_simulation="${output_path}/${protein_root}_${temperature}_${Contacts}.100000"
	input_path="${output_path}/${i}.pdb"
	mv $last_file_from_simulation $input_path
	
	echo "The input path for the next simulation is $input_path"
	awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="NATIVE_FILE"){$2=var}} {print}' cfg_prepare>temp &&  mv temp cfg_prepare  
	#awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' cfg_prepare>temp &&  mv temp cfg_prepare  
	#for diagnostics
	#awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' cfg_prepare
done



