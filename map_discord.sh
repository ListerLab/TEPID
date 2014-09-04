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
<<<<<<< HEAD
            sh process_files.sh $fname $gff &
=======
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
            rm "${fname}.sam"
            rm "{$fname}.bam"
            mv "${fname}_sorted.disc.bed" "${fname}.bed"
            
            echo -e "${blue}Compressing fastq files${NC}"
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq"
            rm "${fname}_1.fastq"
            rm "${fname}_2.fastq"

            echo -e "${blue}Finding TE overlaps${NC}"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $gff > "${fname}_TE_intersections.bed"
            
            echo -e "${green}Finished processing $fname${NC}"
>>>>>>> parent of d558bc7... deleted unnecessary lines
        done
        cd ..
    fi
done
