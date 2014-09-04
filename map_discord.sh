# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -x <path/to/bowtie2/index> -g <path/to/TE_gff>
#   where <proc> is number of processors to use
# Does the following:
#   1. Move into directory
#   2. Maps PE sequencing data using bowtie2
#   3. Converts output to sorted bam file
#   4. Splits discordant reads from bam file into new file
#   5. Deletes concordant reads
#   6. Compresses original fastq files
#   7. Move into next directory and repeat
# Outputs a sorted bedfile containing discordant read alignments and a log file

blue='\033[94m'  # main output
green='\033[92m'  # output for starting / completing files
NC='\033[0m'  # output from samtools etc will be not coloured

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

            echo -e "${blue}Starting mapping${NC}"
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X 3000 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" | samblaster -e -d "${fname}.disc.sam" > /dev/null
            echo -e "${blue}Mapping complete${NC}"
            
            echo -e "${blue}Converting to bam file"
            samtools view -bS "${fname}.disc.sam" | samtools sort -n - "${fname}.sort.disc"

            echo -e "${blue}Converting to bedfile{NC}"
            bedtools bamtobed -bedpe -i "${fname}.sort.disc.bam" > "${fname}.disc.bed"
            
            echo -e "${blue}Sorting bedfile, removing reads mapped to chloroplast or mitochondria${NC}"
            sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${fname}.disc.bed"
            sort -k1,1 -nk2,2 "${fname}.disc.bed" > "${fname}_sorted.disc.bed"
            
            echo -e "${blue}Deleting temp files${NC}"
            rm "${fname}.sort.disc.bam"
            rm "${fname}.disc.bed"
            rm "${fname}.disc.bed.bak"
            rm "${fname}.disc.sam"
            mv "${fname}_sorted.disc.bed" "${fname}.bed"
            
            echo -e "${blue}Compressing fastq files${NC}"
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq"
            rm "${fname}_1.fastq"
            rm "${fname}_2.fastq"

            echo -e "${blue}Finding TE overlaps${NC}"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $gff > "${fname}_TE_intersections.bed"
            
            echo -e "${green}Finished processing $fname${NC}"
        done
        cd ..
    fi
done
