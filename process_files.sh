# Takes sam file containing discordant reads from map_discord.sh and finds TE intersections

# Convert to sorted bam file
samtools view -bS "${1}.disc.sam" | samtools sort -n - "${1}.sort.disc"

# Convert from bam to bedpe
bedtools bamtobed -bedpe -i "${1}.sort.disc.bam" > "${1}.disc.bed"

# Delete reads mapping to chloroplast and mitochondria, rename chromosomes and sort
sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${1}.disc.bed"
sort -k1,1 -nk2,2 "${1}.disc.bed" > "${1}_sorted.disc.bed"

# Remove temp files
rm "${1}.sort.disc.bam"
rm "${1}.disc.bed"
rm "${1}.disc.bed.bak"
rm "${1}.disc.sam"

# Rename bedfile
mv "${1}_sorted.disc.bed" "${1}.bed"

# Find TE intersections
bedtools pairtobed -f 0.1 -type xor -a "${1}.bed" -b $gff > "${1}_TE_intersections.bed"

# Compress fastq files
tar cvfz "${1}.tgz" "${1}_1.fastq" "${1}_2.fastq" --remove-files
