#! /bin/bash

#SBATCH -n 300
#SBATCH -N 6
#SBATCH -J FABG_trunc63_unfolding
#SBATCH -o FABG_trunc63_unfolding.out 
#SBATCH -e FABG_trunc63_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=60000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu

module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 300 ./fold_potential_mpi ./cfg_FABG_trunc63_unfolding > out.txt 32> err.txt 
