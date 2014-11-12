#! /bin/bash
# Created by Tim Stuart

config=  index=  proc=  repo=  yhindex=  size=  genome=  del=  recursive=  helpmsg=  zip=  

# flags that require arguments are followed by :
while getopts C:x:p:c:y:s:g:drhz opt; do
  case $opt in
  C)
      config=$OPTARG
      ;;
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  c)
      repo=${OPTARG%/}
      ;;
  y)
      yhindex=${OPTARG%/}
      ;;
  s)
      size=$OPTARG
      ;;
  g)
      genome=$OPTARG
      ;;
  d)
      del=true
      ;;
  r)
      recursive=true
      ;;
  h)
      helpmsg=true
      ;;
  z)
      zip=true
      ;;
  esac
done
shift $((OPTIND - 1))

if [[ -f $config ]]; then
    . $config
fi

if [ "$helpmsg" == true ]; then
  echo "
  Usage:
    locate.sh [options] -p <proc> -s <size> -x <path/to/bowtie2/index> -y <path/to/yaha/index> -c <path/to/repo> -g <genome>

  Where:
   <proc> is number of processors to use
   <size> is average size of PE fragments sequenced
   <genome> is the organism. Currently supports Arabidopsis (TAIR10 or TAIR9) and Brachypodium.

  Options:
   -C    path to config file. Optional, use instead of command line arguments.
   -d    delete concordantly mapped reads
   -r    run on all subdirectories (recursive)
   -z    gzip input fastq files after mapping
   -h    display help and quit

  Output files:
    * bowtie log file: <name>.log
    * discordant reads bedfile
    * optional concordant reads bamfile
    * split reads bedfile
    * TE insertions bedfile
    * TE deletions bedfile (TE present in Col-0 but not accession)
    * compressed fastq files
    "
  exit
fi

if [ "$zip" == true ]; then
    zip=$zip
else
    zip=false
fi

if [ "$del" == true ]; then
    del=$del
else
    del=false
fi

# Compatable genomes
if [ "$genome" == "TAIR10" ]; then
    gff=$repo/GFF/Arabidopsis/TAIR9_TE.bed
    strip='/Pt/d;/Mt/d;s/chr//g'
elif [ "$genome" == "TAIR9" ]; then
    gff=$repo/GFF/Arabidopsis/TAIR9_TE.bed
    strip='/chrM/d;/chrC/d;s/chr//g'
elif [ "$genome" == "Brachypodium" ]; then
    gff=$repo/GFF/Brachypodium/Brachy_TE_v2.2.bed
    strip='/scaffold_/d;s/Bd//g'
elif [ "$genome" == "hg19" ]; then
    gff=$repo/GFF/Human/hg19.bed
    strip=''
elif [ "$genome" == "hg38" ]; then
    gff=$repo/GFF/Human/hg38.bed
    strip=''
else
    echo "Unsupported genome"
    exit
fi

if [ "$recursive" == true ]; then
  for directory in ./*; do
      if [ -d "$directory" ]; then
          cd $directory
          sh $repo/Code/process_single.sh -p $proc -x $index -c $repo -y $yhindex -s $size -g $genome -d $del -q $strip -z $zip
          cd ..
      fi
  done

else
  sh $repo/Code/process_single.sh -p $proc -x $index -c $repo -y $yhindex -s $size -g $gff -d $del -q $strip -z $zip
fi
