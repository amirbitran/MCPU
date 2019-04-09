#! /bin/bash
#The latest version of the truncate code
#although this doesn't really work with the newest Umbrella platform since stuff is not sent to right path...

#Truncates a desired number of amino acids from an existing protein, and configures everything accordingly

#First argument: number of residues to truncate
#Second argument: Name of protein to truncate (in capital letters)
#Third argument: file root for protein: typically just the protein name in lowercase letters


#Always request an interactive job before running this. Otherwise it may crash
#Also, make sure there is EXACTLY one blank line at the end of the FASTA file
#Modified cfg file is saved as cfg_t$1


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


#linelength=80
#input=104
#if [ "$input" -gt "$linelength" ]; 
#then
#    sed "$ s/\(.\{80\}\)$//" ./mcpu_prep/$2.fasta > ./mcpu_prep/$2_trunc$1.fasta  #first we remove 80
#    remainder=$1-80
#    sed "$ s/\(.\{remainder\}\)$//" ./mcpu_prep/$2_trunc$1.fasta > ./mcpu_prep/$2_trunc$1.fasta  #first we remove 80
#fi




#FIRST STEP WITH FASTA FILE DOESN'T WORK WELL..
#sed "$ s/\(.\{$1\}\)$//" ./mcpu_prep/$2.fasta > ./mcpu_prep/$2_trunc$1.fasta  #removes last N characters from $2.fasta and saves as $2_truncN.fasta



#sed is a text editor
#$ character tells it to work with only the last line
#If we inputted s/.$//, that would remove the last character and replace it with an empty string
#but since we want to remove a variable number of characters, instead of usign ., we use \(.\{$1\}\), which effectively
#tells you to repeat that . N times

#######################################################################
#echo "Requesting interactive session"
#srun -p interact --pty --mem 8000 -t 0-06:00 /bin/bash

#echo "Running C code to generate files for truncated protein"

#UNCOMMENT THE FOLLOWING 3 LINES
cd mcpu_prep
./save_triple $2_trunc$1
cd ..

mkdir ./sim/$2_trunc$1
mkdir ./sim/$2_trunc$1/files
mkdir ./sim/$2_trunc$1/unfolding
mkdir ./sim/$2_trunc$1/replica_tight_spacing

#UNCOMMENT THE FOLLOWING 2 LINES
mv ./mcpu_prep/$2_trunc$1.triple ./sim/$2_trunc$1/files/$2_trunc$1.triple
mv ./mcpu_prep/$2_trunc$1.sctorsion ./sim/$2_trunc$1/files/$2_trunc$1.sctorsion

############################################################################

#Modify secondary structure file
echo "Modifying secondary structure file"
sed "$ s/\(.\{$1\}\)$//" ./sim/$2/files/$2.sec_str > ./sim/$2_trunc$1/files/$2_trunc$1.sec_str  #Remove inputted number of characters from sec_str file


##############################################################

# Truncate PDB file. To do this, start by figuring out how many residues long the protein is
echo "Modifying PDB file"
lastline=$(sed 'x;$!d' <./sim/$2/files/$2.pdb)  #the sed command reads the second to last line of the file (don't want the lat one, since it's "end")



read -a arr <<<$lastline  #Converts last line to an array  

#The fifth element of this array correspodns to the number of amino acids in the protein
#apparently to extract element from array, need to use declare command
declare Size=( "${arr[5]}" )

#Here's what we're gonna do: Read the original pdb file one line at a time. When we reach the line who's amino acid number
#corresponds to Size-$1, we will stop and obtain that line number. This will tell us how many lines we need to keep


#While loop to read line by line...
declare -i res_num=1
declare -i line_num=1

filename="./sim/$2/files/$3_0.100_Emin.pdb"



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
 
(head  -$((${line_num}-1)) ./sim/$2/files/$3_0.100_Emin.pdb; tail -2  ./sim/$2/files/$3_0.100_Emin.pdb)  > ./sim/$2_trunc$1/files/$2_trunc$1.pdb


############################################################################
#Time to modify config files, backbone.c and run.sh

echo "Modifying cfg files"

#First replica

sed "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$3_0.100_Emin.pdb+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.pdb+g" ./src_mpi/cfg_$2 > ./src_mpi/cfg_$2_t$1
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/replica_tight_spacing/$3+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/replica_tight_spacing/$3_t$1+g" ./src_mpi/cfg_$2_t$1
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.triple+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.triple+g" ./src_mpi/cfg_$2_t$1
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.sctorsion+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.sctorsion+g" ./src_mpi/cfg_$2_t$1
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.sec_str+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.sec_str+g" ./src_mpi/cfg_$2_t$1


#Then unfolding
sed "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$3_0.100_Emin.pdb+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.pdb+g" ./src_mpi/cfg_$2_unfolding > ./src_mpi/cfg_$2_t$1_unfolding
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/unfolding/$3+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/unfolding/$3_t$1+g" ./src_mpi/cfg_$2_t$1_unfolding
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.triple+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.triple+g" ./src_mpi/cfg_$2_t$1_unfolding
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.sctorsion+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.sctorsion+g" ./src_mpi/cfg_$2_t$1_unfolding
sed -i "s+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2/files/$2.sec_str+/n/shakfs1/users/amirbitran/ShakhnovichLab/MCPU/sim/$2_trunc$1/files/$2_trunc$1.sec_str+g" ./src_mpi/cfg_$2_t$1_unfolding


#The extension -i tells you to replace in the current file, rather than create a new output file
#In OsX, you need blank quotes after -i. In Odyssey, we do not
#We don't use -i the first time because we haven't created te new file yet



echo "Modifying run.sh"
sed "s+./cfg_$2+./cfg_$2_t$1+g" ./src_mpi/run$2.sh > ./src_mpi/run$2_t$1.sh 
sed "s+./cfg_$2_unfolding+./cfg_$2_t$1_unfolding+g" ./src_mpi/run$2_unfolding.sh > ./src_mpi/run$2_t$1_unfolding.sh 
 


echo "WARNING!!! Please make sure resulting cfg and backbone.c file output simulation to correct directory. Failure to do so may result in data loss...code doesn't consistently do this "
#echo "Exiting interactive session"
#exit 

