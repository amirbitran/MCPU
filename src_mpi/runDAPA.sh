#! /bin/bash

#SBATCH -n 36
#SBATCH -J replicaDAPA
#SBATCH -o replicaDAPA.out 
#SBATCH -e replicaDAPA.err
#SBATCH -p shakhnovich 
#SBATCH --mem=40000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 36 ./fold_potential_mpi ./cfg_DAPA > out.txt 32> err.txt 

