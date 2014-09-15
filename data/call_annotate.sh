path=

while getopts p: opt; do
  case $opt in
  p)
      path=$OPTARG
      ;;
  esac
done
shift $((OPTIND - 1))

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls *intersections.bed);do
        	fname=(${myfile//_TE_intersections.bed/ })
          python "${path}annotate.py" n $fname
          rm "intersections_ordered_${fname}.bed"
        done
        cd ..
    fi
done