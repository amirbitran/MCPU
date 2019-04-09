#! /bin/bash

#SBATCH -n 115
#SBATCH -N 2
#SBATCH -J replica1igd_kbias_0.003
#SBATCH -o replica1igd.out 
#SBATCH -e replica1igd.err
#SBATCH -p shakhnovich 
#SBATCH --mem=50000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 115 ./fold_potential_mpi ./cfg_1igd_3 > out.txt 32> err.txt 
