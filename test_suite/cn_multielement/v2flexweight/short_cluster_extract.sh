#!/bin/bash

# Check if the first argument is "benchmark"
if [[ "$1" == "bm" ]]; then
    echo "Compiling non-GPU accelerated version..."
    g++ -O3 -o extract_clusters extract_clusters.cpp
elif [[ "$1" == "par" ]]; then
    echo "Compiling Parallel-accelerated version..."
    g++ -O3 -o extract_clusters extract_clusters.cpp -fopenmp
else
    echo "Compiling GPU-accelerated version..."
    nvcc -O3 -arch=sm_61 cuda_extract_clusters.cu -o cude_extract_clusters
fi

# Run the compiled program
./cude_extract_clusters setup.in
