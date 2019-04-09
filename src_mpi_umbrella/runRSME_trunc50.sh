#! /bin/bash

#SBATCH -n 225
#SBATCH -N 4
#SBATCH -J replicaRSME_trunc50
#SBATCH -o replicaRSME_trunc50.out 
#SBATCH -e replicaRSME_trunc50.err
#SBATCH -p shakhnovich 
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 225 ./fold_potential_mpi ./cfg_RSME_trunc50 > out.txt 32> err.txt 
