#! /bin/bash

#SBATCH -n 250
#SBATCH -N 4
#SBATCH -J replicaCMK_trunc79
#SBATCH -o replicaCMK_trunc79.out 
#SBATCH -e replicaCMK_trunc79.err
#SBATCH -p shakhnovich  
#SBATCH --mem=60000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 250 ./fold_potential_mpi ./cfg_CMK_trunc79 > out.txt 32> err.txt 
