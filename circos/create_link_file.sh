# usage:
# sh create_link_file.sh -p <path/to/circos/data/dir>

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
        for myfile in $(ls -d insertions*);do
        	fname=(${myfile//insertions_/ })
        	fname=(${fname//.bed/ })
        	sed 1d "insertions_${fname}.bed" > headerless
			awk 'BEGIN {FS=OFS="\t"} {print "chr"$1,$2,$3,"chr"$7,$8,$9}' headerless > "circos_${fname}.txt"
			rm headerless
			cp "circos_${fname}.txt" $path
        done
        cd ..
    fi
done
