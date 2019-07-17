#! /bin/bash

#SBATCH -n 312
#SBATCH -N 5
#SBATCH -J replicaDHFR
#SBATCH -o replicaDHFR.out 
#SBATCH -e replicaDHFR.err
#SBATCH -p shakhnovich
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 312 ./fold_potential_mpi ./cfg_DHFR > out.txt 32> err.txt