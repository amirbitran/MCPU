#! /bin/bash

#SBATCH -n 45
#SBATCH -J ISPA_unfolding
#SBATCH -o ISPA_unfolding.out 
#SBATCH -e ISPA_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=30000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 45 ./fold_potential_mpi ./cfg_ISPA_unfolding > out.txt 32> err.txt 
