#! /bin/bash

#SBATCH -n 324
#SBATCH -N 6
#SBATCH -J replicaISPA
#SBATCH -o replicaISPA.out 
#SBATCH -e replicaISPA.err
#SBATCH -p shakhnovich 
#SBATCH --mem=85000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 324 ./fold_potential_mpi ./cfg_ISPA > out.txt 32> err.txt 
