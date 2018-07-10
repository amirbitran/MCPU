#! /bin/bash

#SBATCH -n 30
#SBATCH -N 1
#SBATCH -J replicaADK_0.650_multistart
#SBATCH -o replicaADK.out 
#SBATCH -e replicaADK.err
#SBATCH -p shakhnovich
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 30 ./fold_potential_mpi ./cfg_ADK > out.txt 32> err.txt
