#! /bin/bash

#SBATCH -n 165
#SBATCH -N 3
#SBATCH -J replica1igd
#SBATCH -o replica1igd.out 
#SBATCH -e replica1igd.err
#SBATCH -p shakhnovich 
#SBATCH --mem=50000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 165 ./fold_potential_mpi2 ./cfg_1igd > out.txt 32> err.txt 
