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
