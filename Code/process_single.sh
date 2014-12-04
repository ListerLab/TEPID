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
bowtie2 --version

for myfile in $(ls -d *_1.fastq);do

    fname=(${myfile//_1.fastq/ })

    echo "Mapping ${fname}"

    bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X $size\
     -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" \
    | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" \
    | samtools view -bS - > "${fname}.bam"

    yaha -t $proc -x $yhindex -q "${fname}.umap.fastq" -L 11 -H 2000 -M 15 -osh stdout \
    | samblaster -s "${fname}.split.sam" > /dev/null

    echo "Converting to bam"

    samtools view -Sb@ $proc "${fname}.split.sam" > "${fname}.split.bam"
    samtools view -Sb@ $proc "${fname}.disc.sam" > "${fname}.disc.bam"

    rm "${fname}.split.sam"
    rm "${fname}.disc.sam"

    if [ "$zip" == true ]; then
      echo "Zipping fastq files"
      gzip *.fastq
    fi

done
