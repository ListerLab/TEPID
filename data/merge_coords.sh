#! /bin/sh

# usage:
# sh merge_coords.sh -f <intersections_file> -a <accession_name> -n <feature_name>

filename=  acc=  name=

while getopts f:a:n: opt; do
  case $opt in
  f)
      filename=$OPTARG
      ;;
  a)
      acc=$OPTARG
      ;;
  n)
     name=$OPTARG
     ;;
  esac
done
shift $((OPTIND - 1))

# create temp directory
mkdir ./temp

# sort by TE name
sort -r -k10,1 $filename > ./temp/sorted_$filename

# split into just TE column (10)
awk 'BEGIN {FS=OFS="\t"} {print $10}' ./temp/sorted_$filename > ./temp/just_tes

# get unique list
uniq ./temp/just_tes > ./temp/uniq_tes

# split each TE into its own file
while read p; do
  grep $p ./$filename > ./temp/temp_$p
done <./temp/uniq_tes

# sort each temp file by chr, start and bedtools merge file. Need current version of bedtools
for myfile in $(ls -d ./temp/temp_*);do
    sort -k1,1 -nk2,2 -o $myfile $myfile
    bedtools merge -c 4,10,6,7,8,12 -o distinct -i $myfile > "${myfile}_merged"
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

