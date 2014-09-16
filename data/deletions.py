# Finds TE insertions that are in Col-0 but not other accessions.
# Run from directory containing all accession subdirectories with mapped data:
# $ python deletions.py g <path/to/TE_gff>
# Note that this script is quite slow.
# Insertions in Col-0 but not PE-sequenced accesions will be evident by discordant mate pairs that:
#  1. map to the same chromosome
#  2. map to different strands
#  3. overlap with a TE
#  4. are sufficiently far apart to accomodate the overlapping TE

import os
from subprocess import call
from sys import argv


def checkArgs(arg1, arg2):
    """
    arg1 is short arg, eg h
    arg2 is long arg, eg host
    """
    args = argv[1:]
    if arg1 in args:
        index = args.index(arg1)+1
        variable = args[index]
        return variable
    elif arg2 in args:
        index = args.index(arg2)+1
        variable = args[index]
        return variable
    else:
        variable = raw_input("\nEnter {arg2}: ".format(arg2=arg2))
        return variable


def find_deletions(bedfile, TEs_dict, accession_name):
    """
    Identifies TE insertions that are in Col-0 but not in other accessions.
    Some bias towards long TEs due to parameters for identifying discordant reads.
    """
    with open(bedfile, 'r') as bed, open('deletions_{a}.bed'.format(a=accession_name), 'w+') as outfile:
        x = 0
        for line in bed:
            line = line.rsplit()
            chr_1 = int(line[0])
            start_1 = int(line[1])
            end_1 = int(line[2])
            chr_2 = int(line[3])
            start_2 = int(line[4])
            end_2 = int(line[5])
            read = line[6]
            strand_1 = line[8]
            strand_2 = line[9]
            if chr_1 == chr_2 and strand_1 != strand_2:
                TE = te_between(chr_1, end_1, start_2, TEs_dict)
                if TE is not False:
                    te_length = abs(TE['start'] - TE['stop'])
                    distance = abs(end_1 - start_2)
                    if te_length <= distance:
                        # write accession names that have the 'deletion'. Later can merge these and change to which have insertion (ie not the 'deletion')
                        outfile.write("{c}\t{s}\t{e}\t{st}\t{n}\t{c}\t{s}\t{e}\t{id}\t{a}\n".format(c=TE['chrom'],
                                                                                                    s=TE['start'],
                                                                                                    e=TE['stop'],
                                                                                                    st=TE['strand'],
                                                                                                    n=TE['name'],
                                                                                                    id='del_' + accession_name + str(x),
                                                                                                    a=accession_name))
                        x += 1
                    else:
                        pass
                else:
                    pass
            else:
                pass


def te_between(chrom, start, stop, TEs_dict):
    """
    Takes set of genomic coordinates (region between two discordant pairs).
    If there is a TE between these coordinates in the Col-0 reference, return dict containing chrom, start, stop, TE name.
    Else return False.
    """
    for key, value in TEs_dict.items():
        if int(value['chrom']) == chrom:
            if overlap(int(value['start']), int(value['stop']), start, stop) is True:
                return {'name': key, 'chrom': int(value['chrom']), 'start': int(value['start']), 'stop': int(value['stop']), 'strand': value['strand']}
            else:
                pass
        else:
            pass
    return False


def overlap(start1, stop1, start2, stop2):
    """
    Takes two pairs of numbers.
    Returns True if sets of coordinates overlap.
    Assumes coordinates are on same chromosome.
    Order of pairs is not important.
    """
    for y in xrange(start2, stop2):
        if start1 <= y <= stop1:
            return True
        else:
            pass

te_file = checkArgs('p', 'path')

with open(te_file, 'r') as infile:
    TEs_dict = {}
    for line in infile:
        line = line.rsplit()
        TEs_dict[line[4]] = {'chrom': line[0], 'start': line[1], 'stop': line[2], 'strand': line[3], 'family': line[4], 'superfamily': line[6]}

for dirs in os.listdir('.'):
    if os.path.isdir(dirs) is True:
        os.chdir(dirs)
        if os.path.isfile('{d}.bed'.format(d=dirs)) is True:
            print dirs
            find_deletions('{d}.bed'.format(d=dirs), TEs_dict, dirs)
            call("sort -k1,1 -nk2,2 deletions_{d}.bed > sorted_deletions_{d}.bed".format(d=dirs), shell=True)
            call("rm deletions_{d}.bed".format(d=dirs))
            os.chdir('..')
        else:
            os.chdir('..')
    else:
        pass
