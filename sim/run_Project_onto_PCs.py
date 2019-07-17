#! /bin/bash

#SBATCH -n 1
#SBATCH -J project_onto_PCs
#SBATCH -o project_onto_PCs.out 
#SBATCH -e project_onto_PCs.err
#SBATCH -p serial_requeue
#SBATCH --mem=8000 
#SBATCH --mail-type=ALL
#SBATCH -t 1-00:00
#SBATCH --mail-user=amirbitran@g.harvard.edu

source new-modules.sh
module load python/3.6.0-fasrc01
python project_onto_PCs.py