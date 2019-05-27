#! /bin/bash

#SBATCH -n 500
#SBATCH -N 8
#SBATCH -J CMK_trunc87_unfolding
#SBATCH -o CMK_trunc87_unfolding.out 
#SBATCH -e CMK_trunc87_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=65000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 500 ./fold_potential_mpi ./cfg_CMK_trunc87_unfolding > out.txt 32> err.txt 
