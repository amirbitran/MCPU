#! /bin/bash

#SBATCH -n 4
#SBATCH -t 2000 
#SBATCH -p general 
#SBATCH --mem=15000 

mpiexec -n 4 ./fold_potential_mpi ./DHFR_replica_test.cfg > out.txt 32> err.txt 
