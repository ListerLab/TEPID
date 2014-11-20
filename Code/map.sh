#! /bin/bash

index=  proc=  yhindex=  size=  recursive=  zip=  repo=  

while getopts x:p:y:s:zc: opt; do
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
  c)
      repo=$OPTARG
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
          sh $repo/Code/process_single.sh -p $proc -x $index -y $yhindex -s $size -z $zip
          cd ..
      fi
  done

else
  sh $repo/Code/process_single.sh -p $proc -x $index -y $yhindex -s $size -z $zip
fi