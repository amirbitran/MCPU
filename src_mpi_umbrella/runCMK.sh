#! /bin/bash

#SBATCH -n 36
#SBATCH -J replicaCMK
#SBATCH -o replicaCMK.out 
#SBATCH -e replicaCMK.err
#SBATCH -p general
#SBATCH --mem=40000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 36 ./fold_potential_mpi ./cfg_CMK > out.txt 32> err.txt 

