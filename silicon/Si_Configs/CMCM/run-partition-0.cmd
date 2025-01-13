#!/bin/bash                                                        
                                                                   
#SBATCH -J block-0                                          
#SBATCH --nodes 1                                                  
#SBATCH --ntasks-per-node 36                                       
#SBATCH -t 1:00:00                                                 
#SBATCH -V                                                         
#SBATCH -o stdoutmsg                                               
#SBATCH -p pbatch                                      
#SBATCH -A iap                                     
time srun -n 36 ./calc_cluster_distance_histograms-mpi training_data training_data s >> run-0.log
