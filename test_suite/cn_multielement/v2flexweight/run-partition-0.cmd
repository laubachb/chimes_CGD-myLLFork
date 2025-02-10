#!/bin/bash                                                        
                                                                   
#SBATCH -J block-0                                          
#SBATCH -N 1                                                       
#SBATCH -n 48                                                      
#SBATCH -t 1:00:00                                                 
#SBATCH -V                                                         
#SBATCH -o stdoutmsg                                               
#SBATCH -p pbatch                                              
#SBATCH -A iap                                             
time srun -n 48 ./calc_cluster_distance_histograms-mpi 000 000 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 001 001 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 002 002 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 003 003 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 004 004 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 005 005 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 006 006 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 007 007 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 008 008 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 009 009 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 010 010 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 011 011 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 012 012 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 013 013 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 014 014 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 015 015 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 016 016 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 017 017 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 018 018 r >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 019 019 r >> run-0.log
