# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -x <path/to/bowtie2/index> -g <path/to/TE_gff> -r <path/to/data/folder>
#   where <proc> is number of processors to use
#   run from directory containing all accession subdirectories
# Does the following:
#   1. Move into subdirectory
#   2. Maps PE sequencing data using bowtie2 local alignment
#   3. Discordant reads are written to sam file, unmapped and soft-clipped reads are written to fastq file (using samblaster)
#   4. Discordant reads are converted into sorted bam file
#   5. Discordant bam file comverted to bedpe file with mate 1 in first position
#   6. Discordant bedpe file is sorted and reads mapping to chloroplast or mitochondrial genomes removed
#   7. Reads where only one mate maps to TE are written to insertions file. Where no reads map to TE, these are written to deletions file
#   8. Python script merges overlapping reads, finds insertion sites
#   9. Python script to create new set of coordinates from deletions file (region between reads), bedtools used to find overlaps with TE (finds deletions)
#   10. Uses yaha to map unaligned or clipped reads using split-read alignment
# X  11. Merges split reads with putative insertion sites to find reads spanning insertion
#   12. Deletes temporary files
#   13. Compresses original fastq files
#   14. Move into next subdirectory and repeat
# Output files:
#  * bowtie metrics file: <name>.log
#  * discordant reads bedfile
#  * split reads bedfile
#  * TE insertions bedfile
#  * TE deletions bedfile (TE present in Col-0 but not accession). Limited to finding TEs >200 bp, <20000 bp
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
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X 250 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" --met-file "${fname}.log" | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" > /dev/null

            echo -e "${blue}Mapping split reads${NC}"
            yaha -t $proc -x $yhindex -q "${fname}.umap.fastq" -L 11 -H 2000 -M 15 -osh "${fname}.split.sam"
            echo -e "${blue}Mapping complete${NC}"            

            echo -e "${blue}Converting to bam file${NC}"
            samtools view -bS "${fname}.disc.sam" | samtools sort -n - "${fname}.sort.disc"
            samtools view -bS "${fname}.split.sam" | samtools sort -n - "${fname}.sort.split"

            echo -e "${blue}Converting to bedfile${NC}"
            bedtools bamtobed -bedpe -mate1 -i "${fname}.sort.disc.bam" > "${fname}.disc.bed"
            bedtools bamtobed -bedpe -mate1 -i "${fname}.sort.split.bam" > "${fname}.split.bed"

            echo -e "${blue}Sorting bedfile, removing reads mapped to chloroplast or mitochondria${NC}"
            sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${fname}.disc.bed"
            sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${fname}.split_unsort.bed"
            sort -k1,1 -nk2,2 "${fname}.disc.bed" > "${fname}.bed"
            sort -k1,1 -nk2,2 "${fname}.split_unsort.bed" > "${fname}.split.bed"

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
            rm "${fname}.split_unsort.bed"

            echo -e "${blue}Compressing fastq files${NC}"
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq" "${fname}.umap.fastq"
            rm "${fname}_1.fastq"
            rm "${fname}_2.fastq"
            rm "${fname}.umap.fastq"

            echo -e "${green}Finished processing $fname${NC}"
        done
        cd ..
    fi
done
