#!/bin/bash
# 02614 - High-Performance Computing, January 2025
# 
# Author: Group 1
#
#BSUB -J poisson_seq_j
#BSUB -o poisson_seq_j_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 10
#BSUB -R "span[hosts=1]"

EXECUTABLE=poisson_j

# Defining parameters
# Number of grid points in each dimension
GRID_SIZES="3 5 10 15 20 25 50 100"
# Starting temp
T_START="-10 -5 0 5 10 15 20 25 30 35"
# Number of iterations
MAX_ITER="1000"
# Tolerance
TOL="1e-3"

for N in $GRID_SIZES
do
    for T in $T_START
    do
        ./$EXECUTABLE $N $MAX_ITER $TOL $T
    done
done
