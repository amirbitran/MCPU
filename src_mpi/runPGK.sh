#! /bin/bash

#SBATCH -n 44
#SBATCH -J replicaPGK
#SBATCH -o replicaPGK.out 
#SBATCH -e replicaPGK.err
#SBATCH -p shakhnovich 
#SBATCH --mem=40000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 44 ./fold_potential_mpi ./cfg_PGK > out.txt 32> err.txt 

