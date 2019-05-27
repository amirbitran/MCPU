#! /bin/bash

#SBATCH -n 75
#SBATCH -N 2
#SBATCH -J replicaMARR_hairpin
#SBATCH -o replicaMARR_hairpin.out 
#SBATCH -e replicaMARR_hairpin.err
#SBATCH -p shakhnovich 
#SBATCH --mem=50000 
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
mpiexec -n 75 ./fold_potential_mpi_orig ./cfg_MARR_hairpin > out.txt 32> err.txt 
