#! /bin/bash

#SBATCH -n 323
#SBATCH -N 6
#SBATCH -J replicaMARR_dimer
#SBATCH -o replicaMARR_dimer.out 
#SBATCH -e replicaMARR_dimer.err
#SBATCH -p shakhnovich 
#SBATCH --mem=95000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 323 ./fold_potential_mpi ./cfg_MARR_dimer > out.txt 32> err.txt 
