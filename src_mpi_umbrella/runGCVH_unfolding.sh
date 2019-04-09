#! /bin/bash

#SBATCH -n 240
#SBATCH -N 4
#SBATCH -J GCVH_unfolding
#SBATCH -o GCVH_unfolding.out 
#SBATCH -e GCVH_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=65000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 240 ./fold_potential_mpi ./cfg_GCVH_unfolding3 > out.txt 32> err.txt 
