#! /bin/bash

#SBATCH -n 36
#SBATCH -J MAP_trunc30_unfolding
#SBATCH -o MAP_trunc30_unfolding.out 
#SBATCH -e MAP_trunc30_unfolding.err
#SBATCH -p general 
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 36 ./fold_potential_mpi ./cfg_MAP_trunc30_unfolding > out.txt 32> err.txt 
