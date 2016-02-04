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
tepid-discover -n <sample_name> -c <mapped> -s <split_mapped> -t <te_bedfile>
```

Where:

  * `<mapped>` is the name of your bam file from bowtie2
  * `<split_mapped` is the name of your split mapped bam file from yaha
  * `<te_bedfile>` is path to the TE annotation bedfile. Annotations included in the repository:  
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

This step is optional and for groups of related samples only, such as a population or generational study. First, a file containing all idenetified variants for the group needs to be generated, with a list of samples that contain each insertion as a column in a bedfile:

```
[chromosome name]	[start position]	[stop position]	[TE name] [comma-separated list of samples containing insertion]
```

For example:
```
chr1	23094	23200	AT1TE69285	Sorbo,Nok-3
```

To do this, the `merge_insertions.py` and `merge_deletions.py` scripts included in the TEpy package, in the `Scripts/` directory, can be used. A list of all sample names is also needed (one sample name on each line of a file). Then, using the merged insertions and deletions files:

```
tepid-refine -t <te_annotation> -i <insertions_file> -d <deletions_file> -a <all_sample_names>
```

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
