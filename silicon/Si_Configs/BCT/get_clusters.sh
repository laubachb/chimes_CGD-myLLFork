#!/bin/bash

# Total frames in INITIAL training set:

nf=1
traj="training_data.xyzf"

g++ -O3 -o extract_clusters extract_clusters.cpp


# This took ~30 min to run

if [ 1 -eq 2 ] ; then

    cp $traj training_data.xyzf
    cp ~/Codes/al_driver-myLLfork/src/helpers.py .
    python -c "import helpers; helpers.break_apart_xyz($nf,\"training_data.xyzf\")"
    rm -f *FORCES*
fi 

if [ 1 -eq 1 ] ; then
    time for i in `ls training_data*xyzf`
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
