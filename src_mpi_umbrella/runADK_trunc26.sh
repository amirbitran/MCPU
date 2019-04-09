#! /bin/bash

#SBATCH -n 350
#SBATCH -N 6
#SBATCH -J replicaADK_trunc26
#SBATCH -o replicaADK_trunc26.out 
#SBATCH -e replicaADK_trunc26.err
#SBATCH -p shakhnovich
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
#. ./PrepareStartingFiles.sh 250 10 25 'ADK_trunc26' 'adk_t26'
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 350 ./fold_potential_mpi ./cfg_ADK_trunc26 > out.txt 32> err.txt
