#! /bin/bash

#SBATCH -n $1
#SBATCH -J $2
#SBATCH -o $2.out 
#SBATCH -e $2.err
#SBATCH -p general 
#SBATCH --mem=15000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n $1 ./fold_potential_mpi ./$2 > out.txt 32> err.txt 
