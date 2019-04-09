#!/bin/bash
#First input is contacts max, second is contacts spacing, third is number of values to sweep over
#Fourth is protein name (ex. 'ADK' or 'ADK_trunc26') that indicates directory where output will be stored
#Fifth is protein root name for output files (ex. 'adk' or 'adk_t26')

#Assumptions:
#The original crystal strucutre filename has the form adk_0.100_Emin.pdb
#There exists a direcotry starting_files (ex. sim/ADK/starting_files) where output will be saved


#an example instance would be . ./PrepareStartingFiles.sh 310 10 31 'ADK' 'adk'



((Contacts_max=$1))
((Contacts_spacing=$2))
((N=$3))
((Contacts_min=Contacts_max-N*Contacts_spacing))
protein_name=$4
protein_root=$5




root='/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/'  #master directory where all PDB files are stored
input_path="$root$protein_name/files/${protein_root}_0.100_Emin.pdb"  #Original (crystal structure) file..assumes of the form ADK/files/adk_0.100_Emin.pdb
output_path="$root$protein_name/starting_files"  #All output (0.pdb, etc) will be saved here


temperature="0.450"


#We will make a new cfg file that will be used to prepare starting files for this particular protein
cfg_file="cfg_prepare_${protein_root}"
cp cfg_prepare $cfg_file


#We use the awk command to edit files...the OFS tells it to use the three tabs as the delimiter
awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="NATIVE_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit native file
awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit structure file
awk -v var="$output_path/$protein_root" '{OFS="\t\t\t"}{if ($1=="PDB_OUT_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit structure file
awk '{OFS="\t\t\t"}{if ($1=="NATIVE_DIRECTORY"){$2="None"}} {print}' $cfg_file>temp && mv temp $cfg_file
awk '{OFS="\t\t\t"}{if ($1=="MC_STEPS"){$2=2000001}} {print}' $cfg_file>temp && mv temp $cfg_file
awk '{OFS="\t\t\t"}{if ($1=="MC_REPLICA_STEPS"){$2=100000000000000000000}} {print}' $cfg_file>temp && mv temp $cfg_file #No replica 
awk '{OFS="\t\t\t"}{if ($1=="MC_PDB_PRINT_STEPS"){$2=2000000}} {print}' $cfg_file>temp && mv temp $cfg_file
awk '{OFS="\t\t\t"}{if ($1=="K_BIAS"){$2=0.5}} {print}' $cfg_file>temp && mv temp $cfg_file
awk '{OFS="\t\t\t"}{if ($1=="NODES_PER_TEMP"){$2=1}} {print}' $cfg_file>temp && mv temp $cfg_file
awk -v var="$temperature" '{OFS="\t\t\t"}{if ($1=="MC_TEMP_MIN"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit structure file


#Had failed to do this previously!
awk -v var="$root$protein_name/files/$protein_name.triple" '{OFS="\t\t\t"}{if ($1=="TRIPLET_ENERGY_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit triplet energy file
awk -v var="$root$protein_name/files/$protein_name.sctorsion" '{OFS="\t\t\t"}{if ($1=="SIDECHAIN_TORSION_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit sidechain torsion file
awk -v var="$root$protein_name/files/$protein_name.sec_str" '{OFS="\t\t\t"}{if ($1=="SECONDARY_STRUCTURE_FILE"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit secondary structure file


for ((i=0; i<N; i++)); do
	#Compute current centerpoint value
	((Contacts=Contacts_max-Contacts_spacing*i))
	
	#Edit cfg file so that simulation uses this centerpoint
	
	awk -v var="$Contacts" '{OFS="\t"}{if ($1=="NUMBER_OF_CONTACTS_MAX"){$2=var}} {print}' $cfg_file>temp && mv temp $cfg_file  #Edit native file
	
	#Run the simulation
	echo "Running simulation at Contacts=$Contacts"
	mpiexec -n 1 fold_potential_mpi $cfg_file
	
	#now, we edit the cfg file so that the next simulation uses the output of the current one to start
	last_file_from_simulation="${output_path}/${protein_root}_${temperature}_${Contacts}.2000000"
	input_path="${output_path}/${i}.pdb"
	mv $last_file_from_simulation $input_path
	
	echo "The input path for the next simulation is $input_path"
	awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="NATIVE_FILE"){$2=var}} {print}' $cfg_file>temp &&  mv temp $cfg_file  
	#awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' $cfg_file>temp &&  mv temp $cfg_file  
	#for diagnostics
	#awk -v var="$input_path" '{OFS="\t\t\t"}{if ($1=="STRUCTURE_FILE"){$2=var}} {print}' $cfg_file
done



