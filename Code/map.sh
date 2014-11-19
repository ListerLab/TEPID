#! /bin/bash

index=  proc=  yhindex=  size=  recursive=  

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
  r)
      recursive=true
      ;;
  esac
done
shift $((OPTIND - 1))

if [ "$recursive" == true ]; then
  for directory in ./*; do
      if [ -d "$directory" ]; then
          cd $directory
          sh process_single.sh -p $proc -x $index -y $yhindex -s $size
          cd ..
      fi
  done

else
  sh process_single.sh -p $proc -x $index -y $yhindex -s $size
fi