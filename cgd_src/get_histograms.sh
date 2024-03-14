#!/bin/bash


system=$1



tag_list=""

# Figure out what frames we have to analyze

for i in `ls training_data_*xyzf`
do
    tag=${i%*.xyzf}
    tag=${tag#*#}

    taglist="${taglist} $tag"
    
done

echo "Compiling the code..."

if [ "$system" = "UM-ARC" ] ; then
    source ~/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod 
elif [ "$system" != "LLNL-LC" ] ; then
    echo "ERROR: System not recognized. Choices are \"UM-ARC\" or \"LLNL-LC"
fi
    
mpiicc -O3 -o calc_cluster_distance_histograms-mpi calc_cluster_distance_histograms-mpi.cpp

# Generate histograms for all the frames
# Going to break the task into 12 blocks, to get everything done in ~1 hr

echo "Building the submit script and submitting..."

jobs_per_block=26 # How many jobs per block
task=0              # Current job
block=-1            # Current block
prev_block=-2       # Previous block


for i in $taglist
do
    if [ `echo "${task} % ${jobs_per_block}" | bc`  == 0 ] ;
    then
    
            let block=block+1
            let prev_block=prev_block+1
    
        
        if [ $prev_block -ge 0 ]; 
        then 
        
            echo "  ... Finished writing block ${prev_block} cmd file ... submitting!"
            echo "  with contents:"
            cat run-partition-${prev_block}.cmd
            sbatch run-partition-${prev_block}.cmd
        fi
        
        # Create the .cmd script
        
        rm -f  run-partition-${block}.cmd        
        echo "#!/bin/bash                                                        "  >> run-partition-${block}.cmd
        echo "                                                                   "  >> run-partition-${block}.cmd
        echo "#SBATCH -J block-${block}                                          "  >> run-partition-${block}.cmd
        echo "#SBATCH --nodes 1                                                  "  >> run-partition-${block}.cmd            
        echo "#SBATCH --ntasks-per-node 36                                       "  >> run-partition-${block}.cmd            
        echo "#SBATCH -t 1:00:00                                                 "  >> run-partition-${block}.cmd        
        echo "#SBATCH -V                                                         "  >> run-partition-${block}.cmd
        echo "#SBATCH -o stdoutmsg                                               "  >> run-partition-${block}.cmd 

        if [ "$system" = "UM-ARC" ] ; then
            echo "#SBATCH -p standard                                            "  >> run-partition-${block}.cmd        
            echo "#SBATCH -A rklinds1                                            "  >> run-partition-${block}.cmd  
            echo "source /home/rklinds/codes/chimes_lsq-myLLfork/modfiles/UM-ARC.mod "  >> run-partition-${block}.cmd          

        elif [ "$system" = "LLNL-LC" ] ; then
                    echo "#SBATCH -p pbatch                                      "  >> run-partition-${block}.cmd        
                    echo "#SBATCH -A pbronze                                     "  >> run-partition-${block}.cmd        
        fi

    fi
    
    echo "time srun -n 36 ./calc_cluster_distance_histograms-mpi $i $i s >> run-${block}.log" >> run-partition-${block}.cmd  
    
    let task=task+1

done

echo "  .... Finished writing the last block ${block} cmd file ... submitting!"
echo "  with contents:"
cat run-partition-${block}.cmd
sbatch run-partition-${block}.cmd; fi

    
