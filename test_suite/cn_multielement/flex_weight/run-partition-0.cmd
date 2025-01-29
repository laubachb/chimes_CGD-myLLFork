#!/bin/bash                                                        
                                                                   
#SBATCH -J block-0                                          
#SBATCH -N 1                                                       
#SBATCH -n 48                                                      
#SBATCH -t 1:00:00                                                 
#SBATCH -V                                                         
#SBATCH -o stdoutmsg                                               
#SBATCH -p pbatch                                              
#SBATCH -A iap                                             
time srun -n 48 ./calc_cluster_distance_histograms-mpi 000 000 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 001 001 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 002 002 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 003 003 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 004 004 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 005 005 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 006 006 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 007 007 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 008 008 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 009 009 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 010 010 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 011 011 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 012 012 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 013 013 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 014 014 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 015 015 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 016 016 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 017 017 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 018 018 s >> run-0.log
time srun -n 48 ./calc_cluster_distance_histograms-mpi 019 019 s >> run-0.log
