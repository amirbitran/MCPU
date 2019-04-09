#! /bin/bash

#SBATCH -n 252
#SBATCH -N 4
#SBATCH -J replicaISPA_trunc87
#SBATCH -o replicaISPA_trunc87.out 
#SBATCH -e replicaISPA_trunc87.err
#SBATCH -p shakhnovich 
#SBATCH --mem=85000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 252 ./fold_potential_mpi ./cfg_ISPA_trunc87 > out.txt 32> err.txt 
