#! /bin/bash

#SBATCH -n 444
#SBATCH -N 7
#SBATCH -J replicaHEMK_trunc15
#SBATCH -o replicaHEMK_trunc15.out 
#SBATCH -e replicaHEMK_trunc15.err
#SBATCH -p shakhnovich  
#SBATCH --mem=70000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 444 ./fold_potential_mpi ./cfg_HEMK_trunc15 > out.txt 32> err.txt 
