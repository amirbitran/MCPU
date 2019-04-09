#! /bin/bash

#SBATCH -n 325
#SBATCH -N 6
#SBATCH -J replicaMARR_trunc32
#SBATCH -o replicaMARR_trunc32.out 
#SBATCH -e replicaMARR_trunc32.err
#SBATCH -p shakhnovich 
#SBATCH --mem=60000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 325 ./fold_potential_mpi ./cfg_MARR_trunc32 > out.txt 32> err.txt 
