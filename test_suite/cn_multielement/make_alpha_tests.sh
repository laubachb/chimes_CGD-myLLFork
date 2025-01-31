#!/bin/bash

# Create the directories and modify the setup.in file
for i in {1..9}; do
    # Create directory alpha_X under tests
    dir="tests/alpha_$i"
    mkdir -p "$dir"
    
    # Copy the required files from flexweight to the new directory
    cp flexweight/*.xyzf "$dir"
    cp flexweight/get_clusters.sh "$dir"
    cp flexweight/extract_clusters.cpp "$dir"
    cp flexweight/calc_cluster_distance_histograms-mpi.cpp "$dir"
    cp flexweight/get_histograms.sh "$dir"
    cp flexweight/helpers.py "$dir"
    cp flexweight/plot_fingerprints.py "$dir"
    cp flexweight/setup.in "$dir"
    
    # Modify the ALPHA value in the setup.in file
    alpha_value=$(echo "scale=4; $i/10" | bc)
    sed -i "s/^ALPHA\s*=.*/ALPHA         = $alpha_value/" "$dir/setup.in"
    
    echo "Directory $dir created and ALPHA set to $alpha_value"
done

# Loop over each alpha_X directory
for dir in tests/alpha_*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        # Submit a job for each directory to run get_clusters.sh
        sbatch <<EOF
#!/bin/bash
#SBATCH -J block-${dir} 
#SBATCH -N 1
#SBATCH -n 12 
#SBATCH -t 1:00:00  
#SBATCH -V                           
#SBATCH -p pdebug 
#SBATCH -A iap

# Change to the directory of the job
cd $dir

# Run the script
./get_clusters.sh
EOF
    fi
done
