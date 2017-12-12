#!/bin/sh
#PBS -l walltime=5:00:00
#PBS -N 2:rep1:m0
#PBS -q gpu
#PBS -l nodes=1:ppn=1
#PBS -m ae
#PBS -M tuf36126@temple.edu
#PBS 

cd $PBS_O_WORKDIR;
mpirun -np 1 ./MIST_0.0.05 -r 2 -n 2 4 -l 1000 -m 0 -q 20 -t 10 -h "(1,2):3" -c 1 1 > log_L_m0; 

wait