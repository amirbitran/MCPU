#! /bin/bash

#SBATCH -n 36
#SBATCH -J NAGA_unfolding
#SBATCH -o NAGA_unfolding.out 
#SBATCH -e NAGA_unfolding.err
#SBATCH -p general 
#SBATCH --mem=50000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 36 ./fold_potential_mpi ./cfg_NAGA_unfolding > out.txt 32> err.txt 
