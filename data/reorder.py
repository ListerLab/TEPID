# usage:
# python reorder.py a <accession> f <feature>
# where <feature> is 'gene' or 'TE'
# Puts TE read in second position, adds column telling which mate is the one mapped to TE (input file has mate 1 in first position)

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


def overlap(start1, stop1, start2, stop2):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    for y in xrange(start2, stop2+1):
        if start1 <= y <= stop1:
            return True
        else:
            pass


def get_len(infile):
    """returns number of lines in file and all lines as part of list"""
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    return i, lines


def reorder(insert_file, reordered_file):
    """
    Reorder columns so that TE read is in second position.
    """
    with open(insert_file, 'r') as infile, open(reordered_file, 'w+') as outfile:
        i, lines = get_len(infile)
        for x in range(i):
            line = lines[x]
            field = line.rsplit()
            read1 = {'chrom': field[0], 'start': field[1], 'stop': field[2], 'strand': field[8]}
            read2 = {'chrom': field[3], 'start': field[4], 'stop': field[5], 'strand': field[9]}
            remaining = [field[10], field[11], field[12], field[13], field[14], field[15], field[6]]  # read name, TE reference coordinates, TE name, TE family
            te_coords = {'chrom': field[10], 'start': field[11], 'stop': field[12], 'strand': field[13], 'name': field[14]}
            if overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                te_read = read1
                dna_read = read2
                mate = 1  # te read is mate 1 in paired data
            elif overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                te_read = read2
                dna_read = read1
                mate = 2
            else:
                raise Exception('check coords')
            outfile.write('{chr1}\t{start1}\t{stop1}\t{strand1}\t{strand2}\t{remain}\t{mate}\n'.format(chr1=dna_read['chrom'],
                                                                                                       start1=dna_read['start'],
                                                                                                       stop1=dna_read['stop'],
                                                                                                       strand1=dna_read['strand'],
                                                                                                       strand2=te_read['strand'],
                                                                                                       remain='\t'.join(remaining),
                                                                                                       mate=mate))

accession = checkArgs('a', 'accession')
f = checkArgs('f', 'feature')

reorder('{acc}_{feat}_intersections.bed'.format(acc=accession, feat=f),
        'intersections_ordered_{feat}_{acc}.bed'.format(acc=accession, feat=f))
