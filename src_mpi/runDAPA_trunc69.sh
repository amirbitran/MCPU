#! /bin/bash

#SBATCH -n 36
#SBATCH -J replicaDAPA_trunc69
#SBATCH -o replicaDAPA_trunc69.out 
#SBATCH -e replicaDAPA_trunc69.err
#SBATCH -p shakhnovich 
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 36 ./fold_potential_mpi ./cfg_DAPA_trunc69 > out.txt 32> err.txt 

