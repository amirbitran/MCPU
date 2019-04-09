#! /bin/bash

#SBATCH -n 125
#SBATCH -N 2
#SBATCH -J replicaMARR_trunc63
#SBATCH -o replicaMARR_trunc63.out 
#SBATCH -e replicaMARR_trunc63.err
#SBATCH -p shakhnovich 
#SBATCH --mem=50000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 125 ./fold_potential_mpi ./cfg_MARR_trunc63 > out.txt 32> err.txt 
