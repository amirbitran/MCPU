#! /bin/bash
#Truncates a desired number of amino acids from MARR, and configures everything accordingly
# You enter the number of amino acids to truncate as first argument
#Always request an interactive job before running this. Otherwise it may crash
#Also, make sure there is EXACTLY one blank line at the end of the FASTA file
#Modified cfg file saved as cfg_t$1 and modified backbone.c saved as backbone_t$1.c

#NOTE; THIS ONLY LETS YOU TRUNCATE UP TO 80 AMINO ACIDS due to how fasta file is formatted
#and the fact that there are 80 characters per line...may want to change this in the future


if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    exit
fi


echo "Please make sure there is EXACTLY one blank line at the end of the FASTA file"

##################################################

#Remove last N characters from {}.fasta. The following works, use it: 

echo "Modifying FASTA file"




sed "$ s/\(.\{$1\}\)$//" ./mcpu_prep/MARR.fasta > ./mcpu_prep/MARR_trunc$1.fasta  #removes last N characters from MARR.fasta and saves as MARR_truncN.fasta



#sed is a text editor
#$ character tells it to work with only the last line
#If we inputted s/.$//, that would remove the last character and replace it with an empty string
#but since we want to remove a variable number of characters, instead of usign ., we use \(.\{$1\}\), which effectively
#tells you to repeat that . N times

#######################################################################
#echo "Requesting interactive session"
#srun -p interact --pty --mem 8000 -t 0-06:00 /bin/bash

#echo "Running C code to generate files for truncated protein"
cd mcpu_prep
./save_triple MARR_trunc$1
cd ..

mkdir ./sim/MARR_trunc$1
mkdir ./sim/MARR_trunc$1/files
mkdir ./sim/MARR_trunc$1/out
mkdir ./sim/MARR_trunc$1/replica_tight_spacing
mkdir ./sim/MARR_trunc$1/unfolding


mv ./mcpu_prep/MARR_trunc$1.triple ./sim/MARR_trunc$1/files/MARR.triple
mv ./mcpu_prep/MARR_trunc$1.sctorsion ./sim/MARR_trunc$1/files/MARR.sctorsion

############################################################################

#Modify secondary structure file
echo "Modifying secondary structure file"
sed "$ s/\(.\{$1\}\)$//" ./sim/MARR/files/MARR.sec_str > ./sim/MARR_trunc$1/files/MARR.sec_str  #Remove inputted number of characters from sec_str file


##############################################################

# Truncate PDB file. To do this, start by figuring out how many residues long the protein is
echo "Modifying PDB file"
lastline=$(sed 'x;$!d' <./sim/MARR/files/MARR.pdb)  #the sed command reads the second to last line of the file (don't want the lat one, since it's "end")



read -a arr <<<$lastline  #Converts last line to an array  

#The fifth element of this array correspodns to the number of amino acids in the protein
#apparently to extract element from array, need to use declare command
declare Size=( "${arr[5]}" )

#Here's what we're gonna do: Read the pdb file one line at a time. When we reach the line who's amino acid number
#corresponds to Size-$1, we will stop and obtain that line number. This will tell us how many lines we need to keep


#While loop to read line by line...
declare -i res_num=1
declare -i line_num=1

filename="./sim/MARR/files/marr_0.100_Emin.pdb"



while read -r line 
do readLine=$line
	if [[ ${line} != "TER" ]] && [[ ${line} != "END" ]];  #For some reason, TER and END were being counted as arrays with residue number 1...This ensures that these lienes are not counted as arrays. Note that this line may not work in all shells
	then 
		read -a arr <<<$line
	fi
	#arr=($line)
	declare res_num=( "${arr[4]}" )
	#echo $res_num
	#declare atom_num=( "${arr[1]}" )
	if [ "$res_num" -le $(( ${Size} - $1 )) ]  
	then
		#echo $"${arr[0]}"
		#echo $line
		line_num=$((line_num+1))   #this will stop adding as soon as the desired condition is reached, plus 1
	fi

done < $filename
#echo $line_num

#The following writes the PDB file onto a new file, keeping only up to line_num-1, AND the last two lines (TER, END)
 
(head  -$((${line_num}-1)) ./sim/MARR/files/marr_0.100_Emin.pdb; tail -2  ./sim/MARR/files/marr_0.100_Emin.pdb)  > ./sim/MARR_trunc$1/files/marr_truncated$1.pdb


############################################################################
#Time to modify config files, backbone.c and run.sh

echo "Modifying cfg file"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/marr_0.100_Emin.pdb+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/marr_truncated$1.pdb+g" ./src_mpi/cfg_MARR > ./src_mpi/cfg_MARR_t$1
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/replica_tight_spacing/marr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/replica_tight_spacing/marr_t$1+g" ./src_mpi/cfg_MARR_t$1
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.triple+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.triple+g" ./src_mpi/cfg_MARR_t$1
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sctorsion+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sctorsion+g" ./src_mpi/cfg_MARR_t$1
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sec_str+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sec_str+g" ./src_mpi/cfg_MARR_t$1


#The extension -i tells you to replace in the current file, rather than create a new output file
#In OsX, you need blank quotes after -i. In Odyssey, we do not
#We don't use -i the first time because we haven't created te new file yet




echo "Modifying run.sh"
sed -i "s+./cfg_MARR+./cfg_MARR_t$1+g" ./src_mpi/runMARR.sh > ./src_mpi/runMARR_t$1.sh
 




############################################################################
#Same thing for the unfolding simulations
echo "Modifying cfg file"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/marr_0.100_Emin.pdb+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/marr_truncated$1.pdb+g" ./src_mpi/cfg_MARR_unfolding > ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/unfolding/marr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/unfolding/marr_t$1+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.triple+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.triple+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sctorsion+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sctorsion+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sec_str+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sec_str+g" ./src_mpi/cfg_MARR_t$1_unfolding


#The extension -i tells you to replace in the current file, rather than create a new output file
#In OsX, you need blank quotes after -i. In Odyssey, we do not
#We don't use -i the first time because we haven't created te new file yet



#The following does not seem to work for some reason?
echo "Modifying run.sh"
sed -i "s+./cfg_MARR_unfolding+./cfg_MARR_t$1_unfolding+g" ./src_mpi/runMARR_unfolding.sh > ./src_mpi/runMARR_t$1_unfolding.sh
 



############################################################################
#Time to get things ready to run :D
echo "Loading modules"
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01



#THE FOLLOWIGN IS SOMETHING YOU NEED TO FIX WHEN YOU EDIT THE CODE
#Also, if you truncate more than 99 amino acids, please change char std_file[100]; to char std_file[101] (or 120, as I have been doing) in backbone
echo "WARNING!!! Please make sure resulting cfg and backbone.c file output simulation to correct directory. Failure to do so may result in data loss...code doesn't consistently do this "
#echo "Exiting interactive session"
#exit 

