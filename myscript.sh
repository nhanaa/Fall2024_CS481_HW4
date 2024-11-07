#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load openmpi/4.1.4-gcc11

mpic++ -g -Wall -O -o hw4 hw4.cpp

mpirun -n 8 hw4 5000 5000 8 /scratch/$USER
mpirun -n 8 hw4 5000 5000 8 /scratch/$USER
mpirun -n 8 hw4 5000 5000 8 /scratch/$USER

mpirun -n 10 hw4 5000 5000 10 /scratch/$USER
mpirun -n 10 hw4 5000 5000 10 /scratch/$USER
mpirun -n 10 hw4 5000 5000 10 /scratch/$USER

mpirun -n 16 hw4 5000 5000 16 /scratch/$USER
mpirun -n 16 hw4 5000 5000 16 /scratch/$USER
mpirun -n 16 hw4 5000 5000 16 /scratch/$USER

mpirun -n 20 hw4 5000 5000 20 /scratch/$USER
mpirun -n 20 hw4 5000 5000 20 /scratch/$USER
mpirun -n 20 hw4 5000 5000 20 /scratch/$USER

module purge
