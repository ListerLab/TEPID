for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d insertions*);do
        	sed 1d "insertions_${directory}.bed" > headerless
			awk 'BEGIN {FS=OFS="\t"} {print "chr"$1,$2,$3,"chr"$7,$8,$9}' "insertions_${directory}.bed" > "circos_${directory}.txt"
			rm headerless
        done
        cd ..
    fi
done
