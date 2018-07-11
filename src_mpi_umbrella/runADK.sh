#! /bin/bash

#SBATCH -n 6
#SBATCH -N 1
#SBATCH -J replicaADK_0.650_multistart
#SBATCH -o replicaADK.out 
#SBATCH -e replicaADK.err
#SBATCH -p shakhnovich
#SBATCH --mem=35000 
#SBATCH --mail-type=ALL
#SBATCH -t 5-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 6 ./fold_potential_mpi ./cfg_ADK > out.txt 32> err.txt
