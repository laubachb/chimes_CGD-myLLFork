#!/bin/bash

#SBATCH -J GPU_cluster
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:00:30
#SBATCH -p spgpu
#SBATCH --gres=gpu:1 
#SBATCH -A eecs570s001w25_class
#SBATCH -V
#SBATCH -o stdoutmsg

module load cuda
rm -rf *.txt
nvcc -O3 -arch=sm_61 cuda_extract_clusters.cu -o cude_extract_clusters
srun -n 1 ./cude_extract_clusters setup.in
 