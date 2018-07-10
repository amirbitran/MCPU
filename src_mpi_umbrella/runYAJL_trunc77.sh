#! /bin/bash

#SBATCH -n 36
#SBATCH -J replicaYAJL_trunc77
#SBATCH -o replicaYAJL_trunc77.out 
#SBATCH -e replicaYAJL_trunc77.err
#SBATCH -p shakhnovich 
#SBATCH --mem=15000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 36 ./fold_potential_mpi ./cfg_YAJL_trunc77 > out.txt 32> err.txt 
