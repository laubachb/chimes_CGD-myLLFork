#!/bin/bash

#SBATCH -J GPU_cluster
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -p spgpu
#SBATCH --gres=gpu:1
#SBATCH -A eecs570s001w25_class
#SBATCH -V
#SBATCH -o stdoutmsg
#SBATCH --mem-per-cpu=5G

module load cuda
rm -rf *.txt
nvcc -O3 -arch=sm_61 test_cuda.cu -o cuda_extract_clusters
time srun -n 1 ./cuda_extract_clusters setup.in
