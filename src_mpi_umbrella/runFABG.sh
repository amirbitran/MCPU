#! /bin/bash

#SBATCH -n 475
#SBATCH -N 8
#SBATCH -J replicaFABG
#SBATCH -o replicaFABG.out 
#SBATCH -e replicaFABG.err
#SBATCH -p shakhnovich 
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 475 ./fold_potential_mpi ./cfg_FABG > out.txt 32> err.txt 
