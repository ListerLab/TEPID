# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -x <path/to/bowtie2/index> -g <path/to/TE_gff> -r <path/to/data/folder>
#   where <proc> is number of processors to use
# Does the following:
#   1. Move into directory
#   2. Maps PE sequencing data using bowtie2
#   3. Converts output to sorted bam file
#   4. Splits discordant reads from bam file into new file
#   5. Annotates TE insertions and deletions from discordant reads
#   6. Deletes temporary files
#   7. Compresses original fastq files
#   8. Move into next directory and repeat
# Output files:
#  * discordant reads bedfile
#  * TE insertions bedfile
#  * TE deletions bedfile (TE present in Col-0 but not accession)
#  * compressed fastq files

blue='\033[94m'  # main output
green='\033[92m'  # output for starting / completing files
NC='\033[0m'  # output from samtools etc will be not coloured

index=  proc=  gff=  pythonpath=

while getopts x:p:g:r: opt; do
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
  r)
      pythonpath=$OPTARG
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
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X 200 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" | samblaster -e -d "${fname}.disc.sam" > /dev/null
            echo -e "${blue}Mapping complete${NC}"

            echo -e "${blue}Converting to bam file${NC}"
            samtools view -bS "${fname}.disc.sam" | samtools sort -n - "${fname}.sort.disc"

            echo -e "${blue}Converting to bedfile${NC}"
            bedtools bamtobed -bedpe -mate1 -i "${fname}.sort.disc.bam" > "${fname}.disc.bed"

            echo -e "${blue}Sorting bedfile, removing reads mapped to chloroplast or mitochondria${NC}"
            sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${fname}.disc.bed"
            sort -k1,1 -nk2,2 "${fname}.disc.bed" > "${fname}.bed"

            echo -e "${blue}Finding TE insertions${NC}"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $gff > "${fname}_TE_intersections.bed"
            python $pythonpath/anotate_ins.py a $fname

            echo -e "${blue}Finding TE deletions${NC}"
            bedtools pairtobed -f 0.1 -type neither -a "${fname}.bed" -b $gff > "${fname}_no_intersections.bed"
            python create_deletion_coords.py b "${fname}_no_intersections.bed" f "${fname}_deletion_coords.bed"
            bedtools intersect -a "${fname}_deletion_coords.bed" -b $gff -wo > "${fname}_deletions_temp.bed"
            python $pythonpath/annotate_del.py a $fname

            echo -e "${blue}Deleting temp files${NC}"
            rm "${fname}.sort.disc.bam"
            rm "${fname}.disc.bed"
            rm "${fname}.disc.bed.bak"
            rm "${fname}.disc.sam"
            rm "${fname}_deletion_coords.bed"
            rm "${fname}_deletions_temp.bed"
            rm "${fname}_no_intersections.bed"
            rm "intersections_ordered_${fname}.bed"
            rm "merged_${fname}.bed"

            echo -e "${blue}Compressing fastq files${NC}"
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq"
            rm "${fname}_1.fastq"
            rm "${fname}_2.fastq"

            echo -e "${green}Finished processing $fname${NC}"
        done
        cd ..
    fi
done
