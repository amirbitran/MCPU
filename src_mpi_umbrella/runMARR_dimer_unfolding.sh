#! /bin/bash

#SBATCH -n 200
#SBATCH -N 4
#SBATCH -J MARR_dimer_unfolding
#SBATCH -o MARR_unfolding.out 
#SBATCH -e MARR_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=85000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 200 ./fold_potential_mpi ./cfg_MARR_dimer_unfolding > out.txt 32> err.txt 
