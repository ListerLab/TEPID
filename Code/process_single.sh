#! /bin/bash

index=  proc=  yhindex=  size=  

while getopts x:p:c:y:s:g:d:q:z: opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  y)
      yhindex=$OPTARG
      ;;
  s)
      size=$OPTARG
      ;;
  esac
done
shift $((OPTIND - 1))

for myfile in $(ls -d *_1.fastq);do

    fname=(${myfile//_1.fastq/ })

    bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X $size\
     -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" --met-file "${fname}.log" \
    | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" \
    | samtools view -bS - > "${fname}.bam"

    yaha -t $proc -x $yhindex -q "${fname}.umap.fastq" -L 11 -H 2000 -M 15 -osh stdout \
    | samblaster -s "${fname}.split.sam" > /dev/null

    samtools view -Sb "${fname}.split.sam" > "${fname}.split.bam"
    samtools view -Sb "${fname}.disc.sam" > "${fname}.disc.bam"

    rm "${fname}.split.sam"
    rm "${fname}.disc.sam"

done
