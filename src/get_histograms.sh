#!/bin/bash

# Setup file parameters

hpc_account=$(grep -oP '^HPC_ACCOUNT\s*=\s*\K.+' setup.in)
system=$(grep -oP '^HPC_SYSTEM\s*=\s*\K.+' setup.in)
hpc_email=$(grep -oP '^HPC_EMAIL\s*=\s*\K.+' setup.in)
email_add=$(grep -oP '^EMAIL_ADD\s*=\s*\K.+' setup.in)
hpc_nodes=$(grep -oP '^HPC_NODES\s*=\s*\K.+' setup.in)
hpc_ppn=$(grep -oP '^HPC_PPN\s*=\s*\K.+' setup.in)
hpc_walltime=$(grep -oP '^HPC_WALLTIME\s*=\s*\K.+' setup.in)
hpc_partition=$(grep -oP '^HPC_PARTITION\s*=\s*\K.+' setup.in)

transformation=$(grep -oP '^TRANSFORMATION\s*=\s*\K.+' setup.in)
jobs_per_block=$(grep -oP '^JOBS_PER_BLOCK\s*=\s*\K.+' setup.in)
src_directory=$(grep -oP '^CGD_SRCDIR\s*=\s*\K.+' setup.in)
cpp_file=${src_directory}/calc_cluster_distance_histograms-mpi.cpp

# Running histogram calculations

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
    
mpiicc -O3 -o calc_cluster_distance_histograms-mpi "$cpp_file"

# Generate histograms for all the frames
# Going to break the task into 12 blocks, to get everything done in ~1 hr

echo "Building the submit script and submitting..."

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
        echo "#SBATCH --nodes ${hpc_nodes}                                       "  >> run-partition-${block}.cmd            
        echo "#SBATCH --ntasks-per-node ${hpc_ppn}                               "  >> run-partition-${block}.cmd            
        echo "#SBATCH -t ${hpc_walltime}                                         "  >> run-partition-${block}.cmd        
        echo "#SBATCH -V                                                         "  >> run-partition-${block}.cmd
        echo "#SBATCH -o stdoutmsg                                               "  >> run-partition-${block}.cmd 

        if [ "$system" = "UM-ARC" ] ; then
            echo "#SBATCH -p ${hpc_partition}                                    "  >> run-partition-${block}.cmd        
            echo "#SBATCH -A ${hpc_account}                                      "  >> run-partition-${block}.cmd          

        elif [ "$system" = "LLNL-LC" ] ; then
                    echo "#SBATCH -p ${hpc_partition}                            "  >> run-partition-${block}.cmd        
                    echo "#SBATCH -A ${hpc_account}                              "  >> run-partition-${block}.cmd        
        fi

    fi
    
    echo "time srun -n ${hpc_ppn} ./calc_cluster_distance_histograms-mpi $i $i $transformation >> cgd_fingerprint.log" >> run-partition-${block}.cmd  
    
    let task=task+1

done

echo "  .... Finished writing the last block ${block} cmd file ... submitting!"
echo "  with contents:"
cat run-partition-${block}.cmd
sbatch run-partition-${block}.cmd

    
