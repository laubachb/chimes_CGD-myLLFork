#!/bin/bash

#SBATCH -J gen_hist
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:01:00
#SBATCH -p debug
#SBATCH -A eecs570s001w25_class
#SBATCH -V
#SBATCH -o stdoutmsg

source ~/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod
mpiicc -O3 -o calc_cluster_distance_histograms-mpi calc_cluster_distance_histograms-mpi.cpp
time srun -n 8 ./calc_cluster_distance_histograms-mpi 000 000 s
 