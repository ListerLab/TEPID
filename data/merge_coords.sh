#! /bin/sh

# usage:
# sh merge_coords.sh -a <accession_name> -n <feature_name>

acc=  name=

while getopts a:n: opt; do
  case $opt in
  a)
      acc=$OPTARG
      ;;
  n)
     name=$OPTARG
     ;;
  esac
done
shift $((OPTIND - 1))

# sort each temp file by chr, start and bedtools merge file. Need current version of bedtools
for myfile in $(ls -d ./temp/temp_*);do
    sort -k1,1 -nk2,2 -o $myfile $myfile
    bedtools merge -c 4,10,6,7,8,9,12,13 -o distinct,distinct,distinct,distinct,distinct,distinct,collapse,collapse -i $myfile > "${myfile}_merged"  # read coords, te ref coords, read name, mate number
done

# cat all the merged TE files. Can be a problem with argument list in ls being too long on some systems, run $ getconf ARG_MAX
for file in $(ls -d ./temp/*_merged);do
  cat $file >> merged_temp
done

# where there are reads on both strands merged (col 4 has +,- or -,+), remove these reads
egrep -v '(+,-)|(-,+)' merged_temp > "merged_${name}_${acc}.bed"

# sort by chr, start
sort -k1,1 -nk2,2 -o "merged_${name}_${acc}.bed" "merged_${name}_${acc}.bed"

# delete temp files
rm -r ./temp
rm merged_temp

