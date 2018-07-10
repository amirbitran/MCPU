#! /bin/bash

#SBATCH -n 10
#SBATCH -J CMK_trunc98_refolding
#SBATCH -o CMK_trunc98_unfolding.out 
#SBATCH -e CMK_trunc98_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 10 ./fold_potential_mpi ./cfg_CMK_trunc98_unfolding > out.txt 32> err.txt 
