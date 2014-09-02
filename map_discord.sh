# Created by Tim Stuart
# Usage:
#   sh map_discord.sh <proc>
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

blue='\033[94m'
NC='\033[0m'

for directory in ./*; do
    if [ -d "$directory" ]; then
        cd $directory
        for myfile in $(ls -d *_1.fastq);do
            fname=(${myfile//_1.fastq/ })
            date
            echo -e "${blue}Starting $fname${NC}"
            echo -e '\tMapping'
            bowtie2 --local -p$1 --fr -q -R5 -N1 -x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9 -X 5000 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" | samblaster -e -d "${fname}.disc.sam" | samtools view -bS - > "${fname}.disc.bam" | tee -a "${fname}.log"
            echo -e '\tFinished mapping $fname'
            echo -e '\tConverting to bedfile'
            bedtools bamtobed -i "${fname}.disc.bam" -bedpe > "${fname}.disc.bed"
            sed -i.bak '/chrC/d;/chrM/d;s/chr//g' "${fname}.disc.bed"
            echo -e '\tSorting bedfile, removing reads mapped to chloroplast or mitochondria'
            sort -k1,1 -nk2,2 "${fname}.disc.bed" > "${fname}_sorted.disc.bed"
            echo -e '\tDeleting temp files'
            rm "${fname}.disc.bed"
            rm "${fname}.disc.bed.bak"
            rm "${fname}.disc.sam"
            rm "${fname}.disc.bam"
            rm "${fname}.sam"
            rm "{$fname}.bam"
            mv "${fname}_sorted.disc.bed" "${fname}.bed"
            echo -e '\tCompressing fastq files'
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq"
            echo -e "${blue}Finished $fname${NC}"
        done
    fi
done
