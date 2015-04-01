*locaTE*
======

Uses paired-end sequencing data to find transposable element insertion points.

Installation
-----

Clone repository

```
git clone git@github.com:timoast/locaTE.git
```

Install requirements

```
pip install -r requirements.txt
```

Run setup.py

```
python setup.py install
```

Usage
-----

Run the mapping script. This is added to your path during the installation.

```
locate_map.sh -x <path/to/bowtie2/index> \
              -p <number_processors> \
              -y <path/to/yaha/index> \
              -s <approximate_fragment_size> \
              -r <recursive> (optional) \
              -z <gzip_fastq_files> (optional)
```

This will give you the following files:

* <name>.bam
* <name>.bam.bai
* <name>.split.bam
* <name>.umap.fastq (this will be compressed if you selected the -z option)

Next go to the directory containing your bam files

```
pylocate -n <sample_name> \
         -c <mapped> \
         -s <split_mapped> \
         -t <te_bedfile>
```

Where:

  * `<mapped>` is the name of your bam file from bowtie2
  * `<split_mapped` is the name of your split mapped bam file from yaha
  * `<te_bedfile>` is path to the TE bedfile included in the repository. Currently supports:  
      - *Arabidopsis thaliana* (TAIR9 and TAIR10)
      - *Brachypodium distachyon*
      - *Homo sapiens* (hg19)

Output files:

  * TE insertions bedfile
  * TE deletions bedfile (TE present in reference but not sample)

---
Required Tools
-------------

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= v2.1.0
* [yaha](https://github.com/GregoryFaust/yaha) >= v0.1.82
* [samtools](http://www.htslib.org/download/) >= v1.1
* [samblaster](https://github.com/GregoryFaust/samblaster) >= v0.1.19
* [bedtools](http://bedtools.readthedocs.org/en/latest/) >= v2.21


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
