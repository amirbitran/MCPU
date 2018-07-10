#! /bin/bash

#SBATCH -n 45
#SBATCH -J replicaP53
#SBATCH -o replicaP53.out 
#SBATCH -e replicaP53.err
#SBATCH -p general 
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 45 ./fold_potential_mpi ./cfg_P53 > out.txt 32> err.txt 
