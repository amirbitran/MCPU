#! /bin/bash

#SBATCH -n 18
#SBATCH -J replicaCRP_unfolding
#SBATCH -o replicaCRP_unfolding.out 
#SBATCH -e replicaCRP_unfolding.err
#SBATCH -p general 
#SBATCH --mem=15000 
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
source new-modules.sh
module load gcc/6.1.0-fasrc01 openmpi/2.0.1-fasrc01
mpiexec -n 18 ./fold_potential_mpi ./cfg_CRP_unfolding > out.txt 32> err.txt 
