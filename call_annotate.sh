for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *intersections.bed);do
        	fname=(${myfile//_TE_intersections.bed/ })
            python annotate.py n $fname
            rm "intersections_ordered_${fname}.bed"
        done
        cd ..
    fi
done