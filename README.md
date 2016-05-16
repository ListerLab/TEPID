*TEPID*
======

Transposable element polymorphism identification

Installation
-----

```
git clone git@github.com:ListerLab/TEPID.git
cd ./TEPID
pip install -r requirements.txt
python setup.py install
```

Usage
-----

Step 1: Mapping
----

```
tepid-map -x <path/to/bowtie2/index> \
         -p <number_processors> \
         -y <path/to/yaha/index> \
         -s <approximate_fragment_size> \
         -r <recursive> (optional) \
         -z <gzip_fastq_files> (optional)
```

This will look for two files named `[name]_1.fastq` and `[name]_2.fastq`, and map these using the number of processors specified in `-p`. These files must be present in the current directory, or in direcories immediately below the current directory if the `-r` option is used.

This will give you the following files:

* [name].bam
* [name].bam.bai
* [name].split.bam
* [name].umap.fastq (this will be compressed if you selected the `-z` option)

The name of these output files will come from the name of the input fastq files

Next go to the directory containing your bam files

Step 2: TE variant discovery
----

```
usage: tepid-discover [-h] [--version] [-k] [-d | -i] [--strict] [--mask MASK]
                      [-D DISCORDANT] [-p PROC] -n NAME -c CONC -s SPLIT -t TE

TEPID -- transposable element polymorphism identification

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -k, --keep            keep all intermediate files
  -d, --deletions       find deletions only
  -i, --insertions      find insertions only
  --strict              Report high-confidence variants only
  --mask MASK           Mask chromosomes in comma separated list or file
  -D DISCORDANT, --discordant DISCORDANT
                        Supply discordant reads bam file
  -p PROC, --proc PROC  number of processors
  -n NAME, --name NAME  sample name
  -c CONC, --conc CONC  bam file from bowtie2
  -s SPLIT, --split SPLIT
                        split reads bam file from yaha
  -t TE, --te TE        TE annotation bedfile
```

The following TE annotations for use with TEPID are included in the repository:  
  - *Arabidopsis thaliana* (TAIR9 and TAIR10)
  - *Brachypodium distachyon*
  - *Homo sapiens* (hg19)

Output files:

  * TE insertions bedfile
  * TE deletions bedfile (TE present in reference but not sample)
  * File containing names of reads providing evidence for insertions
  * File containing names of reads providing evidence for deletions

Step 3: Refinement and genotyping
----

```
usage: tepid-refine [-h] [--version] [-k] [-i INSERTIONS] [-d DELETIONS]
                    [-p PROC] -t TE -n NAME -c CONC -s SPLIT -a ALL_SAMPLES

TEPID -- refine TE insertion and deletion calls

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -k, --keep            keep all intermediate files
  -i INSERTIONS, --insertions INSERTIONS
                        File containing collapsed TE insertions for all
                        samples in population
  -d DELETIONS, --deletions DELETIONS
                        File containing collapsed TE deletions for all samples
                        in population
  -p PROC, --proc PROC  number of processors
  -t TE, --te TE        TE annotation bedfile
  -n NAME, --name NAME  sample name
  -c CONC, --conc CONC  bam file from bowtie2
  -s SPLIT, --split SPLIT
                        split reads bam file from yaha
  -a ALL_SAMPLES, --all_samples ALL_SAMPLES
                        List of all sample names
```

This step is optional and for groups of related samples only, such as a population or generational study. First, a file containing all idenetified variants for the group needs to be generated, with a list of samples that contain each insertion as a column in a bedfile:

```
[chromosome name]	[start position]	[stop position]	[TE name] [comma-separated list of samples containing insertion]
```

For example:
```
chr1	23094	23200	AT1TE69285	Sorbo,Nok-3
```

To do this, the `merge_insertions.py` and `merge_deletions.py` scripts included in the TEpy package, in the `Scripts/` directory, can be used. A list of all sample names is also needed (one sample name on each line of a file).

---
Required Tools
-------------

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= v2.1.0
* [yaha](https://github.com/GregoryFaust/yaha) >= v0.1.82
* [samtools](http://www.htslib.org/download/) >= v1.1
* [samblaster](https://github.com/GregoryFaust/samblaster) >= v0.1.19
* [bedtools](http://bedtools.readthedocs.org/en/latest/) >= v2.25.0


**Python requirements**

* [python](https://www.python.org) v2.7
* [numpy](http://www.numpy.org/)
* [pybedtools](http://pythonhosted.org/pybedtools/)
* [pysam](http://pysam.readthedocs.org/en/latest/)

---
License
-------

This software is licensed under the GNU General Public License. See ./LICENSE
for further details. If this license is incompatible with your open source
project, talk to us (raise an issue) and we will negotiate suitable re- or
dual-licensing.
