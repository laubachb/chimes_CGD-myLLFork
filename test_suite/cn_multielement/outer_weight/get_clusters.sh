#!/bin/bash

# Total frames in INITIAL training set:
rm -f *.txt *.hist run*
#nf=20
#traj="/g/g10/laubach2/chimes_CGD-myLLFork/test_suite/cn_multielement/test/data.xyzf"

# Read N_FRAMES and TRAJ_PATH from setup.in
nf=$(grep -oP '^N_FRAMES\s*=\s*\K\d+' setup.in)
traj=$(grep -oP '^TRAJ_PATH\s*=\s*\K.+' setup.in)

# Total frames in INITIAL training set:
echo "Using N_FRAMES: $nf"
echo "Using TRAJ_PATH: $traj"

g++ -O3 -o extract_clusters extract_clusters.cpp

#for i in `ls training_data_#019.xyzf`
# This took ~30 min to run

if [ 1 -eq 2 ] ; then

    cp $traj training_data.xyzf
    cp /g/g10/laubach2/codes/al_driver-myLLfork/src/helpers.py .
    python -c "import helpers; helpers.break_apart_xyz($nf,\"training_data.xyzf\")"
    rm -f *FORCES*
fi 

if [ 1 -eq 1 ] ; then
    time for i in `ls training_data_*xyzf`
    #for i in `ls training_data_#019.xyzf`
    do
        tag=${i%*.xyzf}
        tag=${tag#*#}
    
        echo "Working on tag: $tag"
    
        cp $i test.xyz
    
        time ./extract_clusters setup.in
        
        # Rename _clu-*.txt files with the tag
        for j in 2b_clu-*.txt 3b_clu-*.txt 4b_clu-*.txt
        do
            mv "$j" "${tag}.${j}"
        done

        # Process and rename clutypes files
        for j in clutypes_*b.txt
        do
            if [[ $j == clutypes_2* ]]; then
                sed -i '1d' "$j"  # Remove first line
            elif [[ $j == clutypes_3* ]]; then
                sed -i '1,3d' "$j"  # Remove first three lines
            elif [[ $j == clutypes_4* ]]; then
                sed -i '1,6d' "$j"  # Remove first six lines
            fi
            mv "$j" "${tag}.${j}"
        done

        # Now paste the renamed clutypes files with the corresponding clu files
        for j in ${tag}.2b_clu-*.txt ${tag}.3b_clu-*.txt ${tag}.4b_clu-*.txt
        do
            # Extract 2b, 3b, or 4b from filename
            clu_type=$(echo "$j" | grep -o '[234]b')

            # Correctly find the corresponding clutypes file
            clutype_file="${tag}.clutypes_${clu_type}.txt"

            if [[ -f "$clutype_file" ]]; then
                paste "$j" "$clutype_file" > temp.txt && mv temp.txt "$j"
            else
                echo "Warning: Matching file $clutype_file not found for $j"
            fi
        done
        
    done
fi

# # Run histogram script at the end
# echo "Running get_histograms.sh..."
# time ./get_histograms.sh
