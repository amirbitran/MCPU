#! /bin/bash

#SBATCH -n 1
#SBATCH -J PrepareStartingFiles
#SBATCH -o PrepareStartingFiles.out 
#SBATCH -e PrepareStartingFiles.err
#SBATCH -p shakhnovich
#SBATCH --mem=25000 
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu
module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02
. ./PrepareStartingFiles.sh 250 10 25 'ADK_trunc26' 'adk_t26'
. ./PrepareStartingFiles.sh 330 10 26 'FABG_trunc63' 'fabg_t63'


