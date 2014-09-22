Population epigenetics
======================

Schmitz et al. analyzed genetic and epigenetic variation in over 150 different *Arabidopsis* accessions using paired-end Illumina sequencing and whole-genome bisulfite sequencing.  

They found several thousand differentially methylated regions between the accessions. PCR analysis of a small, randomly selected subset of these differentially methylated regions revealed that approximately 15% are coupled to genetic changes such as transposable element insertions.

However, this finding was not pursued further. I aim to use the published paired-end sequencing data to identify structural changes between the accessions, and look for correlations between these changes and the previously identified differentially methylated regions.

This repository will provide all the information needed to repeat my analysis, including:  
1. Downloading data from the SRA.  
2. Splitting `.sra` files into paired `.fastq` files.  
3. Mapping of the paired `.fastq` files to the *Arabidopsis* reference genome.  
4. Extracting discordant reads from the alignment files.  
5. Using the discordant reads to identify potential transposable element insertion sites.  
6. Correlating the putative transposable element insertion sites with the published set of differentially methylated regions from Schmitz et al. 2013 *Nature*.  

### Full reference:  
Schmitz, R. J., Schultz, M. D., Urich, M. A., Nery, J. R., Pelizzola, M., Libiger, O., et al. (2013). Patterns of population epigenomic diversity. Nature, 495(7440), 193–198. doi:10.1038/nature11968


### Required tools:
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://samtools.sourceforge.net)
* [samblaster](https://github.com/GregoryFaust/samblaster)
* [bedtools](http://bedtools.readthedocs.org/en/latest/) v2.21 or newer
