#! /bin/bash

index=  proc=  yhindex=  size=  recursive=  zip=  

while getopts x:p:y:s:z opt; do
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
  z)
      zip=true
      ;;
  esac
done
shift $((OPTIND - 1))

if [ "$zip" == true ]; then
  zip=$zip
else
  zip=false
fi

if [ "$recursive" == true ]; then
  for directory in ./*; do
      if [ -d "$directory" ]; then
          cd $directory
          sh process_single.sh -p $proc -x $index -y $yhindex -s $size -z $zip
          cd ..
      fi
  done

else
  sh process_single.sh -p $proc -x $index -y $yhindex -s $size -z $zip
fi