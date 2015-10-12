#! /usr/bin/env python

import pysam
import pybedtools
import os
from tepy import tepy


"""

Need to supply insertions file and deletions file separately, as deletions lines are inverse in poly file

With condensed te file for all accessions, examine sites in accessions without the identified TE indel for:
  -  concordant reads spanning site (evidence of no indel)
  -  if no concordant reads spanning, look for split reads and return split_names
  -  extract reads from fastq file and remap, intersect with TEs and add to list of accessions if there is evidence of TE insertions (any number of reads)

"""


def readNames(names):
    n = []
    with open(names, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            n.append(line[0])
    return n


def getOtherLines(acc, infile):
    """
    Return list of accessions for which TE variant was not identified, with name of TE(s) and coordinates
    """
    regions = []
    with open(infile, 'r') as inf:
        for line in inf:
            line = line.rsplit()
            accessions = line[-1].split(',')
            inverse_acc = [i for i in acc if i not in accessions]
            regions.append((line, inverse_acc))
    return regions


def find_reads(coords, bam):
    r = []
    bam_reads = bam.fetch(coords[0], coords[1], coords[2])
    for i in bam_reads:
        r.append(i.qname)
    bam.reset()
    if len(r) > 0:
        return r
    else:
        return False


def extract_reads(bam, name_indexed, names, acc):
    header = bam.header.copy()
    out_name = 'extracted_reads_{}.bam'.format(acc)
    out_bam = pysam.Samfile(out_name, 'wb', header=header)
    for name in names:
        iterator = name_indexed.find(name)
        for x in iterator:
            out_bam.write(x)
    out_bam.close()
    return out_name


def check_te_overlaps(te, bamfile, te_list):
    intersections = pybedtools.BedTool(bamfile).bam_to_bed().intersect(te, wb=True)
    reads = []
    for r in intersections:
        if r[-3] in te_list:
            reads.append(r[3])
        else:
            pass
    if len(reads) > 0:
        return reads
    else:
        return False


def get_last_id(acc, indel):
    with open("{i}_reads_{a}.txt".format(i=indel, a=acc), 'r') as f:
        for line in f:
            pass
    return int(line.rsplit()[0][1:])


def write_te(te_file_name, read_file_name, data, read_names, iterator):
    with open(te_file_name, 'a+') as te_file, open(read_file_name, 'a+') as read_file:
        coords = [str(x) for x in i[0]]
        te_file.write("\t".join(coords)+"\t"+i[2]+"\n")
        read_file.write(">"+iterator+"\t"+",".join(read_names)+"\n")


def process_missed(data, indel, concordant, split_alignments, name_indexed, acc):
    for i in data:
        coords = (i[0][0], int(i[0][1]), int(i[0][2]))
        te_list = i[0][-2].split(",")
        if acc in i[1]:
            # if find_reads(coords, concordant) is False:  # this will be a problem, maybe just look for split reads
            split_names = find_reads(coords, split_alignments)
            if split_names is not False:
                extracted = extract_reads(split_alignments, name_indexed, split_names, acc)
                read_names = check_te_overlaps(te, extracted, te_list)
                if read_names not False:
                    read_file_name = "second_pass_reads_{t}_{a}.txt".format(t=indel, a=acc)
                    te_file_name = "second_pass_deletions_{t}_{a}.bed".format(t=indel, a=acc)
                    try:
                        iterator
                    except NameError:
                        iterator = get_last_id(acc, indel)
                    else:
                        iterator += 1
                    write_te(te_file_name, read_file_name, i[0], read_names, iterator)
                else:
                    pass
            else:
                pass
            # else:
            #     pass
        else:
            pass


def refine(options):
    te = pybedtools.BedTool(options.te).sort()
    names = readNames(options.all_samples)
    insertions = getOtherLines(names, options.insertions)
    deletions = getOtherLines(names, options.deletions)  # format ([data], [inverse_accessions])
    for acc in os.listdir('.'):
        if os.path.isdir(acc) is True:
            os.chdir(acc)
            if os.path.isfile('deletions_{}.bed'.format(acc)) is True:
                print "Processing "+acc

                conc = acc+"_filtered.bam"
                split = acc+".split.bam"
                tepy.check_bam(conc, options.proc)
                tepy.check_bam(split, options.proc)
                concordant = pysam.AlignmentFile(conc, 'rb')
                split_alignments = pysam.AlignmentFile(split, 'rb')
                name_indexed = pysam.IndexedReads(split_alignments)
                name_indexed.build()

                print "  deletions"
                process_missed(deletions, "deletions", concordant, split_alignments, name_indexed, acc)
                print "  insertions"
                process_missed(insertions, "insertions", concordant, split_alignments, name_indexed, acc)

                os.chdir('..')
            else:
                os.chdir('..')
        else:
            pass


if __name__ == "__main__":
    from argparse import ArgumentParser
    import pkg_resources

    version = pkg_resources.require("TEpy")[0].version

    parser = ArgumentParser(description='Perform second pass over population data to catch false negatives')
    parser.add_argument('-k', '--keep', help='keep all intermediate files', action='store_true', required=False, default=False)
    parser.add_argument('-p', '--proc', help='number of processors', required=False, default=1, type=int)
    parser.add_argument('-t', '--te', help='TE annotation bedfile', required=True)
    parser.add_argument('--version', action='version', version='%(prog)s '+str(version))
    parser.add_argument('-i', '--insertions', help='File containing collapsed TE insertions for all samples in population', required=True)
    parser.add_argument('-d', '--deletions', help='File containing collapsed TE deletions for all samples in population', required=True)
    parser.add_argument('-a', '--all_samples', help='List of all sample names', required=True)
    options = parser.parse_args()
    refine(options)
