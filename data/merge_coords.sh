#! /bin/sh

filename=  acc=

while getopts f:a: opt; do
  case $opt in
  f)
      filename=$OPTARG
      ;;
  a)
      acc=$OPTARG
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

# split each TE into its own file (for each strand)
while read p; do
  grep $p ./$filename > ./temp/temp_$p
done <./temp/uniq_tes

# sort each temp file by chr, start and bedtools merge file
for myfile in $(ls -d ./temp/temp_*);do
    sort -k1,1 -nk2,2 -o $myfile $myfile
    bedtools merge -c 4,10,6,7,8,12 -o distinct -i $myfile> "${myfile}_merged"  # specify stramd in bedtools
done

# cat all the merged TE files
cat ./temp/*_merged > merged_$acc.bed

# sort by chr, start
sort -k1,1 -nk2,2 -o merged_$acc.bed merged_$acc.bed

# delete temp files
rm -r ./temp

