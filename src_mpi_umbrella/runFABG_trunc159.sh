#! /bin/bash

#SBATCH -n 175
#SBATCH -N 3
#SBATCH -J replicaFABG_trunc159
#SBATCH -o replicaFABG_trunc159.out 
#SBATCH -e replicaFABG_trunc159.err
#SBATCH -p shakhnovich 
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 175 ./fold_potential_mpi ./cfg_FABG_trunc159 > out.txt 32> err.txt 
