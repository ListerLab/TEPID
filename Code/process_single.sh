#! /bin/sh

blue='\033[94m'  # main output
green='\033[92m'  # output for starting / completing files
NC='\033[0m'  # output from samtools etc will be not coloured

index=  proc=  repo=  yhindex=  size=  gff=  keep=  strip=  zip=

while getopts x:p:c:y:s:g:k:q:z: opt; do
  case $opt in
  x)
      index=$OPTARG
      ;;
  p)
      proc=$OPTARG
      ;;
  c)
      repo=$OPTARG
      ;;
  y)
      yhindex=$OPTARG
      ;;
  s)
      size=$OPTARG
      ;;
  g)
      gff=$OPTARG
      ;;
  k)
      keep=$OPTARG
      ;;
  q)
      strip=$OPTARG
      ;;
  z)
      zip=$OPTARG
      ;;
  esac
done
shift $((OPTIND - 1))

for myfile in $(ls -d *_1.fastq);do
    fname=(${myfile//_1.fastq/ })
    date

    echo -e "${green}Processing $fname${NC}"

    echo -e "${blue}Starting mapping${NC}"
    bowtie2 --local --dovetail -p$proc --fr -q -R5 -N1 -x $index -X $size -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" --met-file "${fname}.log" | samblaster -e -d "${fname}.disc.sam" -u "${fname}.umap.fastq" > $keep

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
    sed -i.bak $strip "${fname}.disc.bed"
    sed -i.bak $strip "${fname}.split_unsort.bed"
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
    if [ "$keep" != "/dev/null" ]; then
      read mean std <<< $(head -5000 $keep | python $repo/Code/calc_mean.py)
    else
      mean=False  std=False
    fi
    python $repo/Code/create_deletion_coords.py b "${fname}_no_TE_intersections.bed" f "${fname}_deletion_coords.bed" m $mean d $std s $size
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

    if [ "$zip" == true ]; then
      echo -e "${blue}Compressing fastq files${NC}"
      gzip "${fname}_1.fastq" "${fname}_2.fastq" "${fname}.umap.fastq"
    fi

    echo -e "${green}Finished processing ${fname}${NC}"
done
