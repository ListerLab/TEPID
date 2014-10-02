#! /bin/sh

# Created by Tim Stuart
# Usage:
#   sh map_discord.sh -p <proc> -x <path/to/bowtie2/index> -y <path/to/yaha/index> -r <path/to/repo>
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
#  * bowtie log file: <name>.log
#  * discordant reads bedfile
#  * split reads bedfile
#  * TE insertions bedfile
#  * TE deletions bedfile (TE present in Col-0 but not accession). Limited to finding TEs >200 bp, <20000 bp
#  * compressed fastq files

# To do:
#  * merge split read and discordant read data

blue='\033[94m'  # main output
green='\033[92m'  # output for starting / completing files
NC='\033[0m'  # output from samtools etc will be not coloured

index=  proc=  repo=  yhindex=

while getopts x:p:g:r:y: opt; do
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
            # if we want to integrate this with other analysis, should have option to keep concordantly mapped files
            bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X 250 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" --met-file "${fname}.log" | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" > /dev/null

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
            sort -k1,1 -nk2,2 "${fname}.split_unsort.bed" > "${fname}.split.bed"  # need to do something with split reads

            echo -e "${blue}Finding insertions${NC}"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $repo/GFF/TAIR9_TE.bed > "${fname}_TE_intersections.bed"
            bedtools pairtobed -f 0.1 -type xor -a "${fname}.bed" -b $repo/GFF/TAIR10_genes.bed > "${fname}_gene_intersections.bed"
            python $repo/data/reorder.py a $fname f TE
            sh $repo/data/merge_coords.sh -f "intersections_ordered_TE_${fname}.bed" -a $fname -n TE
            python $repo/data/annotate_ins.py a $fname f TE
            python $repo/data/reorder.py a $fname f gene
            sh $repo/data/merge_coords.sh -f "intersections_ordered_gene_${fname}.bed" -a $fname -n gene
            python $repo/data/annotate_ins.py a $fname f gene

            echo -e "${blue}Finding deletions${NC}"
            bedtools pairtobed -f 0.1 -type neither -a "${fname}.bed" -b $gff/TAIR9_TE.bed  > "${fname}_no_te_intersections.bed"
            bedtools pairtobed -f 0.1 -type neither -a "${fname}.bed" -b $gff/TAIR10_genes.bed> "${fname}_no_gene_intersections.bed"
            python $repo/data/create_deletion_coords.py b "${fname}_no_te_intersections.bed" f "${fname}_deletion_coords.bed"
            bedtools intersect -a "${fname}_deletion_coords.bed" -b $gff -wo > "${fname}_deletions_temp.bed"
            python $repo/data/annotate_del.py a $fname

            echo -e "${blue}Deleting temp files${NC}"
            # need to fix some of these filenames
            rm "${fname}.sort.disc.bam"
            rm "${fname}.disc.bed"
            rm "${fname}.disc.bed.bak"
            rm "${fname}.disc.sam"
            rm "${fname}_deletion_coords.bed"
            rm "${fname}_deletions_temp.bed"
            rm "${fname}_no_te_intersections.bed"
            rm "${fname}_no_gene_intersections.bed"
            rm "intersections_ordered_te_${fname}.bed"
            rm "intersections_ordered_gene_${fname}.bed"
            rm "merged_gene_${fname}.bed"
            rm "merged_te_${fname}.bed"
            rm "${fname}.split_unsort.bed"
            rm "${fname}.split_unsort.bed.bak"

            echo -e "${blue}Compressing fastq files${NC}"
            gzip "${fname}_1.fastq" "${fname}_2.fastq" "${fname}.umap.fastq"

            echo -e "${green}Finished processing $fname${NC}"
        done
        cd ..
    fi
done
