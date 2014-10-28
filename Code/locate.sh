#! /bin/sh
# Created by Tim Stuart

index=  proc=  repo=  yhindex=  size=  genome=  keep=  helpmsg=  

# required flags are followed by :
while getopts x:p:c:y:s:g:kh opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  c)
      repo=$OPTARG
      ;;
  y)
      yhindex=$OPTARG
      ;;
  s)
      size=$OPTARG
      ;;
  g)
      genome=$OPTARG
      ;;
  k)
      keep=$OPTARG
      ;;
  r)
      recursive=$OPTARG
      ;;
  h)
      helpmsg=true
      ;;
  esac
done
shift $((OPTIND - 1))

if [ "$helpmsg" == true ]; then
  echo "\tUsage:

  \tlocate.sh [options] -p <proc> -s <size> -x <path/to/bowtie2/index> -y <path/to/yaha/index> -c <path/to/repo> -g <genome>

  \twhere:
 \t\t<proc> is number of processors to use
 \t\t<size> is average size of PE fragments sequenced
 \t\t<genome> is the organism. Currently supports Arabidopsis and Brachypodium.

 \tOptions:
 \t\t-k <path>    keep concordantly mapped reads, store at <path>
 \t\t-r           run on all subdirectories (recursive)

\tOutput files:
\t\t* bowtie log file: <name>.log
\t\t* discordant reads bedfile
\t\t* split reads bedfile
\t\t* TE insertions bedfile
\t\t* TE deletions bedfile (TE present in Col-0 but not accession). Limited to finding TEs >200 bp, <20000 bp
\t\t* compressed fastq files"
  exit
fi

if [ "$genome" == "Arabidopsis" ]; then
    gff=$repo/GFF/Arabidopsis/TAIR9_TE.bed
elif [ "$genome" == "Brachypodium" ]; then
    gff=$repo/GFF/Brachypodium/Brachy_TE_v2.2.bed
else
    echo "Unsupported genome"
    exit
fi

if [ "$keep" ]; then
    keep=$keep
else
    keep=/dev/null
fi

if [ "$recursive" ]; then
  for directory in ./*; do
      if [ -d "$directory" ]; then
          cd $directory
          sh $repo/Code/process_single.sh -p $proc -x $index -c $repo -y $yhindex -s $size -g $genome -k $keep
          cd ..
      fi
  done

else
  sh $repo/Code/process_single.sh -p $proc -x $index -c $repo -y $yhindex -s $size -g $gff -k $keep
fi
