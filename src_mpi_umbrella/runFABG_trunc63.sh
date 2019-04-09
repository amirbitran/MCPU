#! /bin/bash

#SBATCH -n 325
#SBATCH -N 6
#SBATCH -J replicaFABG_trunc63
#SBATCH -o replicaFABG_trunc63.out 
#SBATCH -e replicaFABG_trunc63.err
#SBATCH -p shakhnovich 
#SBATCH --mem=60000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu

module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 325 ./fold_potential_mpi ./cfg_FABG_trunc63 > out.txt 32> err.txt 
