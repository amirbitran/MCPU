#! /bin/bash

#SBATCH -n 18
#SBATCH -J ADK_unfolding
#SBATCH -o replicaCMK.out 
#SBATCH -e replicaCMK.err
#SBATCH -p general 
#SBATCH --mem=40000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 18 ./fold_potential_mpi ./cfg_ADK_unfolding > out.txt 32> err.txt 

