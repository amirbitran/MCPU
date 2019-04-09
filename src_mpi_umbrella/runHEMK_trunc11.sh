#! /bin/bash

#SBATCH -n 518
#SBATCH -N 9
#SBATCH -J replicaHEMK_trunc11
#SBATCH -o replicaHEMK_trunc11.out 
#SBATCH -e replicaHEMK_trunc11.err
#SBATCH -p shakhnovich  
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 518 ./fold_potential_mpi ./cfg_HEMK_trunc11 > out.txt 32> err.txt 
