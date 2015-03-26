#! /bin/bash

index=  proc=  yhindex=  size=  

while getopts x:p:y:s:z: opt; do
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
  z)
      zip=$OPTARG
      ;;
  esac
done
shift $((OPTIND - 1))

samtools --version
samblaster --version

for myfile in $(ls -d *_1.fastq);do

    fname=(${myfile//_1.fastq/ })

    echo "Mapping ${fname}"
    bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X $size\
     -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" \
    | samblaster -e -u "${fname}.umap.fastq" \
    | samtools view -b - > "${fname}.bam"

    yaha -t $proc -x $yhindex -q "${fname}.umap.fastq" -L 11 -H 2000 -M 15 -osh stdout \
    | samblaster -s "${fname}.split.sam" > /dev/null

    echo "Converting to bam"
    samtools view -Sb@ $proc "${fname}.split.sam" > "${fname}.split.bam"
    rm "${fname}.split.sam"

    echo "Sorting alignment"
    samtools sort -@ $proc "${fname}.bam" "${fname}_sorted"
    rm "${fname}.bam"

    echo "Indexing alignment"
    samtools index "${fname}_sorted.bam"

    if [ "$zip" == true ]; then
      echo "Zipping fastq files"
      gzip *.fastq
    fi

done
