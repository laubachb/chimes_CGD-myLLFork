#!/bin/bash

# Total frames in INITIAL training set:

nf=$(grep -oP '^NFRAMES\s*=\s*\K\d+' setup.in)
traj=$(grep -oP '^TRAJPATH\s*=\s*\K.+' setup.in)
src_directory=$(grep -oP '^CGD_SRCDIR\s*=\s*\K.+' setup.in)
work_directory=$(grep -oP '^WORKING_DIR\s*=\s*\K.+' setup.in)
cpp_file=${src_directory}/extract_clusters.cpp

g++ -O3 -o extract_clusters "$cpp_file"

# This took ~30 min to run

if [ 1 -eq 1 ] ; then

    cp $traj ${work_directory}/training_data.xyzf
    cp "${src_directory}/helpers.py" .
    python -c "import helpers; helpers.break_apart_xyz($nf,\"training_data.xyzf\")"
    rm -f *FORCES*
fi 

if [ 1 -eq 1 ] ; then
    time for i in `ls training_data_*xyzf`
    do

        tag=${i%*.xyzf}
        tag=${tag#*#}
    
        echo "Working on tag: $tag"
    
        cp $i test.xyz
    
        time ./extract_clusters
        
        for j in `ls 2b_clu-*.txt 3b_clu-*.txt 4b_clu-*.txt`
        do
            mv $j ${tag}.${j}
        done

    done
fi
