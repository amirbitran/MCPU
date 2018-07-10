#! /bin/bash
#Truncates a desired number of amino acids from DHFR, and configures everything accordingly
# You pass N, the number of amino acids to truncate, as first input
#Always request interactive job before running this
#Modified cfg file saved as cfg_t$1 and modified backbone.c saved as backbone_t$1.c







##################################################

#Remove last N characters from 4DFR.fasta. The following works, use it: 

echo "Modifying FASTA file"
sed "$ s/\(.\{$1\}\)$//" ./mcpu_prep/4DFR.fasta > ./mcpu_prep/4DFR_trunc$1.fasta  #removes last N characters from 4DFR.fasta and saves as $4DFR_truncN.fasta
#sed is a text editor
#$ character tells it to work with only the last line
#If we inputted s/.$//, that would remove the last character and replace it with an empty string
#but since we want to remove a variable number of characters, instead of usign ., we use \(.\{$1\}\), which effectively
#tells you to repeat that . N times

#######################################################################
#echo "Requesting interactive session"
#srun -p interact --pty --mem 8000 -t 0-06:00 /bin/bash

echo "Running C code to generate files for truncated protein"
./mcpu_prep/save_triple 4DFR

mkdir ./sim/DHFR_trunc$1
mkdir ./sim/DHFR_trunc$1/files
mkdir ./sim/DHFR_trunc$1/out
mkdir ./sim/DHFR_trunc$1/replica



mv ./mcpu_prep/4DFR.triple ./sim/DHFR_trunc$1/files
mv ./mcpu_prep/4DFR.sctorsion ./sim/DHFR_trunc$1/files

############################################################################

#Modify secondary structure file
echo "Modifying secondary structure file"
sed "$ s/\(.\{$1\}\)$//" ./sim/DHFR/files/4DFR.sec_str > ./sim/DHFR_trunc$1/files/4DFR.sec_str  #Remove inputted number of characters from sec_str file


##############################################################

# Truncate PDB file. To do this, start by figuring out how many residues long the protein is
echo "Modifying PDB file"
lastline=$(sed 'x;$!d' <./sim/DHFR/files/4DFR.pdb)  #the sed command reads the second to last line of the file (don't want the lat one, since it's "end")



read -a arr <<<$lastline  #Converts last line to an array  

#The fifth element of this array correspodns to the number of amino acids in the protein
#apparently to extract element from array, need to use declare command
declare Size=( "${arr[5]}" )

#Here's what we're gonna do: Read the pdb file one line at a time. When we reach the line who's amino acid number
#corresponds to Size-$1, we will stop and obtain that line number. This will tell us how many lines we need to keep


#While loop to read line by line...
declare -i res_num=1
declare -i line_num=1

filename="./sim/DHFR/files/dhfr_0.100_Emin.pdb"



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
 
(head  -$((${line_num}-1)) ./sim/DHFR/files/dhfr_0.100_Emin.pdb; tail -2  ./sim/DHFR/files/dhfr_0.100_Emin.pdb)  > ./sim/DHFR_trunc$1/files/dhfr_truncated$1.pdb


############################################################################
#Time to modify config files, backbone.c and run.sh

echo "Modifying cfg file"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/files/dhfr_0.100_Emin.pdb+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/files/dhfr_truncated$1.pdb+g" ./src_mpi/cfg > ./src_mpi/cfg_t$1
sed -i "" "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/replica/dhfr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/replica/dhfr_t$1+g" ./src_mpi/cfg_t$1
sed -i "" "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/files/4DFR.triple+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/files/4DFR.triple+g" ./src_mpi/cfg_t$1
sed -i "" "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/files/4DFR.sctorsion+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/files/4DFR.sctorsion+g" ./src_mpi/cfg_t$1
sed -i "" "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/files/4DFR.sec_str+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/files/4DFR.sec_str+g" ./src_mpi/cfg_t$1


#The extension -i tells you to replace in the current file, rather than create a new output file
#In OsX, you need blank quotes after -i . In Odyssey, I believe you do not
#We don't use -i the first time because we haven't created te new file yet


echo "Modifying backbone.c"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/replica/dhfr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/replica/dhfr_t$1+g" ./src_mpi/backbone.c > ./src_mpi/backbone_t$1.c



echo "Modifying run.sh"
sed "s+./fold_potential_mpi+./fold_potential_trunc$1+g" ./src_mpi/run.sh > ./src_mpi/run_t$1.sh
sed -i "" "s+./cfg+./cfg_t$1+g" ./src_mpi/run.sh 



############################################################################
#Compile C Code, to get things ready to run :D
echo "Loading modules"
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01

echo "Compiling C Code"
mpicc -O3 -o fold_potential_trunc$1 backbone_t$1.c -lm 

echo "Exiting interactive session"
exit 

