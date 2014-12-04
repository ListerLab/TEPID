*locaTE*
======

Uses paired-end sequencing data to find transposable element insertion points.

Installation
-----

Clone repository

```
git clone git@github.com:timoast/locaTE.git
```

Usage
-----

Go to the directory containing your bam files

```
python find_te.py -c <mapped> -d <disc_mapped> -s <split_mapped> -t <te_bedfile>
```

Where:

  * `<mapped>` is the name of your mapped bam file
  * `<disc_mapped>` is the name of your discordant mapped bam file
  * `<split_mapped` is the name of your split mapped bam file
  * `<te_bedfile>` is path to the TE bedfile included in the repository. Currently supports:  
      - *Arabidopsis thaliana* (TAIR9 and TAIR10)
      - *Brachypodium distachyon*
      - *Homo sapiens* (hg19)

Output files:

  * bedfiles for split and discordant reads
  * TE insertions bedfile
  * TE deletions bedfile (TE present in reference but not sample)

---
Required Tools
-------------

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= v2.1.0
* [yaha](http://faculty.virginia.edu/irahall/yaha/) >= v0.1.82
* [samtools](http://samtools.sourceforge.net) >= v1.1
* [samblaster](https://github.com/GregoryFaust/samblaster) >= v0.1.19
* [bedtools](http://bedtools.readthedocs.org/en/latest/) >= v2.21


**Python requirements**

* [python](https://www.python.org) v2.7
* [numpy](http://www.numpy.org/)
* [pybedtools](http://pythonhosted.org/pybedtools/)

---
License
-------

This software is licensed under the GNU General Public License. See ./LICENSE
for further details. If this license is incompatible with your open source
project, talk to us (raise an issue) and we will negotiate suitable re- or
dual-licensing.
