#! /bin/bash

#SBATCH -n 200
#SBATCH -N 4
#SBATCH -J replicaCMK_trunc87
#SBATCH -o replicaCMK_trunc87.out 
#SBATCH -e replicaCMK_trunc87.err
#SBATCH -p shakhnovich  
#SBATCH --mem=60000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 200 ./fold_potential_mpi ./cfg_CMK_trunc87 > out.txt 32> err.txt 
