#! /bin/bash

#SBATCH -n 36
#SBATCH -J YAJL_unfolding
#SBATCH -o YAJL_unfolding.out 
#SBATCH -e YAJL_unfolding.err
#SBATCH -p shakhnovich 
#SBATCH --mem=15000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 27 ./fold_potential_mpi ./cfg_YAJL_unfolding > out.txt 32> err.txt 
