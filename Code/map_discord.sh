#! /bin/sh

# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -s <size> -x <path/to/bowtie2/index> -y <path/to/yaha/index> -r <path/to/repo> -g <genome>
#   where:
#  <proc> is number of processors to use
#  <size> is average size of PE fragments sequenced
#  <genome> is the organism. Currently supports Arabidopsis and Brachypodium.

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
#   11. Merges split reads with putative insertion sites to find reads spanning insertion
#   12. Deletes temporary files
#   13. Compresses original fastq files
#   14. Move into next subdirectory and repeat

# Output files:
#  * bowtie log file: <name>.log
#  * discordant reads bedfile
#  * split reads bedfile
#  * TE insertions bedfile
#  * TE deletions bedfile (TE present in Col-0 but not accession). Limited to finding TEs >200 bp, <20000 bp
#  * compressed fastq files

blue='\033[94m'  # main output
green='\033[92m'  # output for starting / completing files
NC='\033[0m'  # output from samtools etc will be not coloured

index=  proc=  repo=  yhindex=  size=  genome=

while getopts x:p:r:y:s:g: opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  r)
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
  esac
done
shift $((OPTIND - 1))

if [ "$genome" == "Arabidopsis" ]; then
    gff = $repo/GFF/Arabidopsis/TAIR9_TE.bed
elif [ "$genome" == "Brachypodium" ]; then
    gff = $repo/GFF/Brachypodium/Brachy_TE_v2.2.bed
else
    echo "Unsupported genome"
    exit
fi

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *_1.fastq);do
            fname=(${myfile//_1.fastq/ })
            date

            echo -e "${green}Processing $fname${NC}"

            echo -e "${blue}Starting mapping${NC}"
            # if we want to integrate this with other analysis, should have option to keep concordantly mapped files
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X $size -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" --met-file "${fname}.log" | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" > /dev/null

            echo -e "${blue}Mapping split reads${NC}"
            yaha -t $proc -x $yhindex -q "${fname}.umap.fastq" -L 11 -H 2000 -M 15 -osh stdout | samblaster -s "${fname}.split.sam" > /dev/null
            echo -e "${blue}Mapping complete${NC}"

            echo -e "${blue}Converting to bam file${NC}"
            samtools view -bS "${fname}.disc.sam" | samtools sort -n - "${fname}.sort.disc"
            samtools view -bS "${fname}.split.sam" | samtools sort -n - "${fname}.sort.split"

            echo -e "${blue}Converting to bedfile${NC}"
            bedtools bamtobed -bedpe -mate1 -i "${fname}.sort.disc.bam" > "${fname}.disc.bed"
            bedtools bamtobed -i "${fname}.sort.split.bam" > "${fname}.split_unsort.bed"

            echo -e "${blue}Sorting bedfile, removing reads mapped to chloroplast or mitochondria${NC}"
            sed -i.bak '/Pt/d;/Mt/d;s/chr//g' "${fname}.disc.bed"  # tair10 annotation, tair9 uses chrC and chrM for chloroplast, mitochondrial genomes. Would need to be changed for other organisms
            sed -i.bak '/Pt/d;/Mt/d;s/chr//g' "${fname}.split_unsort.bed"
            sort -k1,1 -nk2,2 "${fname}.disc.bed" > "${fname}.bed"
            sort -k1,1 -nk2,2 "${fname}.split_unsort.bed" > "${fname}.split.bed"

            echo -e "${blue}Finding insertions${NC}"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $gff > "${fname}_TE_intersections.bed"
            python $repo/Code/reorder.py a $fname f TE
            sort -k10 "intersections_ordered_TE_${fname}.bed" > "intersections_ordered_TE_${fname}_sort.bed"
            mkdir ./temp
            cd ./temp
            python $repo/Code/split_bed_by_gene.py ../"intersections_ordered_TE_${fname}_sort.bed" 9 temp
            cd ..
            sh $repo/Code/merge_coords.sh -a $fname -n TE
            python $repo/Code/annotate_ins.py a $fname f TE
            bedtools merge -c 2,3 -o count_distinct,count_distinct -i "${fname}.split.bed" > "${fname}_merged.split.bed"
            python $repo/Code/filter_split.py $fname
            bedtools intersect -c -a "insertions_TE_${fname}_temp.bed" -b "${fname}_filtered.split.bed" > "insertions_TE_${fname}_split_reads.bed"
            python $repo/Code/separate_breakpoints.py $fname
            bedtools intersect -a single_break_temp.bed -b "${fname}_filtered.split.bed" -wo > single_break.bed
            bedtools intersect -a double_break_temp.bed -b "${fname}_filtered.split.bed" -wo > double_break.bed
            python $repo/Code/annotate_breakpoints.py $fname
            python $repo/Code/separate_reads.py $fname
            sort -k1,1 -nk2,2 "insertions_${fname}_unsorted.bed" > "insertions_${fname}.bed"

            echo -e "${blue}Finding deletions${NC}"
            bedtools pairtobed -f 0.1 -type neither -a "${fname}.bed" -b $gff  > "${fname}_no_TE_intersections.bed"
            python $repo/Code/create_deletion_coords.py b "${fname}_no_TE_intersections.bed" f "${fname}_deletion_coords.bed"
            bedtools intersect -a "${fname}_deletion_coords.bed" -b $gff -wo > "${fname}_deletions_temp.bed"
            python $repo/Code/annotate_del.py a $fname

            echo -e "${blue}Deleting temp files${NC}"
            rm "${fname}.sort.disc.bam"
            rm "${fname}.disc.bed"
            rm "${fname}.disc.bed.bak"
            rm "${fname}.disc.sam"
            rm "${fname}_deletion_coords.bed"
            rm "${fname}_deletions_temp.bed"
            rm "${fname}_no_TE_intersections.bed"
            rm "intersections_ordered_TE_${fname}.bed"
            rm "intersections_ordered_TE_${fname}_sort.bed"
            rm "${fname}.split_unsort.bed"
            rm "${fname}.split_unsort.bed.bak"
            rm "insertions_TE_${fname}_split_reads.bed"
            rm "insertions_TE_${fname}_temp.bed"
            rm "${fname}_merged.split.bed"
            rm "${fname}_filtered.split.bed"
            rm single_break_temp.bed
            rm double_break_temp.bed
            rm single_break.bed
            rm double_break.bed
            rm "insertions_${fname}_unsorted.bed"
            rm "insertions_${fname}_temp.bed"
            rm "${fname}.sort.split.bam"
            rm "${fname}.split.sam"
            rm "merged_TE_${fname}.bed"
            rm "${fname}_TE_intersections.bed"

            echo -e "${blue}Compressing fastq files${NC}"
            gzip "${fname}_1.fastq" "${fname}_2.fastq" "${fname}.umap.fastq"

            echo -e "${green}Finished processing $fname${NC}"
        done
        cd ..
    fi
done
