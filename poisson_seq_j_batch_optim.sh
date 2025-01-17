#!/bin/bash
# 02614 - High-Performance Computing, January 2025
# 
# Author: Group 1
#
#BSUB -J poisson_seq_j_opt
#BSUB -o poisson_seq_j_opt_%J.out
#BSUB -q hpcintro
#BSUB -n 1
#BSUB -R "rusage[mem=2048]"
#BSUB -W 10
#BSUB -R "span[hosts=1]"

EXECUTABLE=poisson_j

# Defining parameters
# Number of grid points in each dimension
N="100"
# Starting temp
T="0"
# Number of iterations
MAX_ITER="1000"
# Tolerance
TOL="1e-3"

# Compiler options
# COMPILER_OPTIONS=("O0" "O1" "O2" "O3" "Ofast")
OPT="Ofast"

./$EXECUTABLE $N $MAX_ITER $TOL $T
echo "N=$N, T=$T, OPT=$OPT"
# done