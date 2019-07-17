#! /bin/bash

#SBATCH -n 275
#SBATCH	-N 5
#SBATCH -J replicaGCVH_trunc21
#SBATCH -o replicaGCVH_trunc21.out 
#SBATCH -e replicaGCVH_trunc21.err
#SBATCH -p shakhnovich  
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 275 ./fold_potential_mpi ./cfg_GCVH_trunc21 > out.txt 32> err.txt 