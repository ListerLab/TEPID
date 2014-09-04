# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -x <path/to/bowtie2/index> -g <path/to/TE_gff>
#   where <proc> is number of processors to use
# Does the following:
#   1. Move into directory
#   2. Maps PE sequencing data using bowtie2
#   3. Extracts discordant reads
#   4. Calls process_files.sh and starts mapping next accession

green='\033[92m'
NC='\033[0m'

index=  proc=  gff=

while getopts x:p:g: opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  g)
      gff=$OPTARG
      ;;
  esac
done
shift $((OPTIND - 1))

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *_1.fastq);do
            fname=(${myfile//_1.fastq/ })
            date
            echo -e "${green}Processing $fname${NC}"
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X 3000 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" | samblaster -e -d "${fname}.disc.sam" > /dev/null
            sh process_files.sh $fname &
        done
        cd ..
    fi
done
