#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load openmpi/4.1.4-gcc11

mpic++ -g -Wall -O -o hw4 hw4.cpp

mpirun -n 2 hw4 100 100 2 /scratch/$USER
mpirun -n 4 hw4 100 100 8 /scratch/$USER
mpirun -n 8 hw4 100 100 8 /scratch/$USER
mpirun -n 16 hw4 100 100 16 /scratch/$USER

diff /scratch/$USER/output100.100.2.0 /scratch/$USER/output100.100.4.0
diff /scratch/$USER/output100.100.2.0 /scratch/$USER/output100.100.8.0
diff /scratch/$USER/output100.100.4.0 /scratch/$USER/output100.100.8.0
diff /scratch/$USER/output100.100.8.0 /scratch/$USER/output100.100.16.0

module purge
