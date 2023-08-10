#!/bin/bash

tag_list=""

# Figure out what frames we have to analyze

for i in `ls training_data_*xyzf`
do
    tag=${i%*.xyzf}
    tag=${tag#*#}

    taglist="${taglist} $tag"
    
done

tagarr=($taglist)

echo "Compiling the code..."


g++ -O3 -o CluUQ_similarity CluUQ_similarity.cpp

echo "Running the tasks"

ntags=${#tagarr[@]}

rm -f dist_matrix.dat

time for (( i=0; i<$ntags; i++ ))
do
	itag=`printf "%04d" "$i"`

	echo "Running loop over i= $i"
	for (( j=i+1; j<ntags; j++ ))
	do
		jtag=`printf "%04d" "$j"`
	
		del_2b=`./CluUQ_similarity ${itag}-${itag}.2b_clu-s.hist ${jtag}-${jtag}.2b_clu-s.hist`
		del_3b=`./CluUQ_similarity ${itag}-${itag}.3b_clu-s.hist ${jtag}-${jtag}.3b_clu-s.hist`
		del_4b=`./CluUQ_similarity ${itag}-${itag}.4b_clu-s.hist ${jtag}-${jtag}.4b_clu-s.hist`
		
		del_all=`echo "$del_2b + $del_3b + $del_4b" | bc -l`
		
		echo "$i $j $del_2b $del_3b $del_4b $del_all" >> dist_matrix.dat
	done
	
	echo "" >> dist_matrix.dat
done

    
