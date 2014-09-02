# Created by Tim Stuart
# Usage:
#   sh map_discord.sh <proc>
#   where <proc> is number of processors to use
# Requirements:
#   samblaster: https://github.com/GregoryFaust/samblaster
#   bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#   samtools: http://samtools.sourceforge.net/samtools.shtml
# Does the following:
#   1. Move into directory
#   2. Maps PE sequencing data using bowtie2
#   3. Converts output to sorted bam file
#   4. Splits discordant reads from bam file into new file
#   5. Deletes concordant reads
#   6. Compresses original fastq files
#   7. Move into next directory and repeat

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
            bowtie2 --local -p$1 --fr -q -R5 -N1 -x /dd_stage/userdata/lister/data/genomes/bowtie2_indexes/tair9 -X 2000 -1 "${fname}_1.fastq" -2 "${fname}_2.fastq" -S "${fname}.sam" | tee -a "${fname}.log"
            echo -e '\tFinished mapping $fname'
            echo -e '\tConverting output to sorted bam file'
            samtools view -bS $fname$samname | samtools sort -n - $fname
            echo -e '\tSplitting discordant reads from bam file'
            samtools view -h "${fname}.bam" | samblaster -a -e -d "${fname}.disc.bam"
            rm "${fname}.sam"
            rm "{$fname}.bam"
            echo -e '\tCompressing fastq files'
            tar cvfz "${fname}.tgz" "${fname}_1.fastq" "${fname}_2.fastq"
            echo -e "${blue}Finished $fname${NC}"
        done
    fi
done
