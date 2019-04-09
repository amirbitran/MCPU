#! /bin/bash

#SBATCH -n 200
#SBATCH -N 4
#SBATCH -J replicaFABG_trunc116
#SBATCH -o replicaFABG_trunc116.out 
#SBATCH -e replicaFABG_trunc116.err
#SBATCH -p shakhnovich 
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 200 ./fold_potential_mpi ./cfg_FABG_trunc116 > out.txt 32> err.txt 
