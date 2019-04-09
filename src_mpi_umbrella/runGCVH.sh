#! /bin/bash

#SBATCH -n 300
#SBATCH	-N 5
#SBATCH -J replicaGCVH
#SBATCH -o replicaGCVH.out 
#SBATCH -e replicaGCVH.err
#SBATCH -p shakhnovich  
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 300 ./fold_potential_mpi ./cfg_GCVH > out.txt 32> err.txt 
