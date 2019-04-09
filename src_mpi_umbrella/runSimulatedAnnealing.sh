#! /bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -J Simulated_annealing
#SBATCH -o Simulated_annealing.out 
#SBATCH -e Simulated_annealing.err
#SBATCH -p shakhnovich 
#SBATCH --mem=10000 
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
. ./Simulated_annealing.sh 1igd 1igd