from __future__ import division
import os
from sys import argv
import numpy as np
import pybedtools


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


def filter_lines(feature, max_dist):
    """
    Filters discordant reads and removes
    reads mapped to mitochondria and chloroplast
    use in pybedtools.filter()
    """
    start1 = int(feature[2])
    start2 = int(feature[4])
    chr1 = feature[0]
    chr2 = feature[3]
    unwanted_chroms = ['Mt', 'Pt', 'chrM', 'chrC']
    if chr1 in unwanted_chroms or chr2 in unwanted_chroms:
        return False
    elif abs(start1 - start2) > max_dist or chr1 != chr2:
        return True
    else:
        return False


def filter_lines_split(feature):
    """
    Removes reads mapped to mito and chl
    use in pybedtools.filter()
    """
    unwanted_chroms = ['Mt', 'Pt', 'chrM', 'chrC']
    if feature[0] in unwanted_chroms:
        return False
    else:
        return True


def remove_chr(feature):
    """
    use in pybedtools.each()
    """
    feature[0] = str(feature[0]).strip('chrBd')
    feature[3] = str(feature[3]).strip('chrBd')
    return feature


def _create_te_dict(infile):
    """
    Create dictionary where key is TE name
    value is another dictionary where key is 
    read number (arbitrary), value is line in original
    bedfile
    """
    TE_dict = {}
    for line in infile:
        fields = line.rsplit()
        line = line.strip().split('\t')
        name = fields[9]
        try:
            TE_dict[name]
        except KeyError:
            TE_dict[name] = {1: line}
        else:
            x = len(TE_dict[name])
            TE_dict[name][x+1] = line
    return TE_dict


def merge_te_coords(infile, outfile):
    """
    takes file containing reads coordinates
    that overlap annotated TEs and creates a 
    new file with merged coordinates
    """
    with open(infile, 'r') as inp:
        TE_dict = _create_te_dict(inp)
    with open(outfile, 'a+') as outf:
        for name in TE_dict.keys():
            _modify_coords(TE_dict[name], outf)
            for key, value in TE_dict[name].items():
                chrom = value[0]
                start = value[1]
                stop = value[2]
                strand = value[5]
                ref = '\t'.join(value[3])
                reads = ','.join(value[4])
                mates = ','.join(value[6])
                if len(strand) == 1:
                    outf.write('{ch}\t{sta}\t{sto}\t{str}\t{name}\t{ref}\t{reads}\t{mates}\n'.format(ch=chrom,
                                                                                                     sta=start,
                                                                                                     sto=stop,
                                                                                                     str=strand[0],
                                                                                                     name=name,
                                                                                                     ref=ref,
                                                                                                     reads=reads,
                                                                                                     mates=mates))
                else:
                    pass


def _modify_coords(inp, outf):
    """
    Go through each read mapped to TE and merge
    if they overlap, modify dictionary with new
    information
    """
    if len(inp) > 1:
        skips = []  # keys that have already been merged
        for x in inp.keys():
            if x not in skips:
                line = inp[x]
                chrom1 = line[0]
                start1 = int(line[1])
                stop1 = int(line[2])
                reads = [line[11]]
                strands = [line[3]]
                ref_coords = line[5:9]
                mates = [line[12]]
                _merge(chrom1, start1, stop1, inp, x, strands, reads, mates, ref_coords, outf, skips)
            else:
                pass
    else:
        _no_merging(1, inp)


def _no_merging(x, inp):
    """
    modify format of unmerged read coordinates
    """
    line = inp[x]
    chrom1 = line[0]
    start1 = int(line[1])
    stop1 = int(line[2])
    reads = [line[11]]
    strands = [line[3]]
    ref_coords = line[5:9]
    mates = [line[12]]
    inp[x] = [chrom1, start1, stop1, ref_coords, reads, strands, mates]


def _merge(chrom1, start1, stop1, d, x, strands, reads, mates, ref_coords, outf, skips):
    """
    merge overlapping read coordinates where TE name is the same
    """
    for key in d.keys():
        if key == x or key in skips:
            pass
        else:
            line = d[key]
            chrom2 = line[0]
            start2 = int(line[1])
            stop2 = int(line[2])
            read2 = line[11]
            strand2 = line[3]
            mate2 = line[12]
            if chrom1 == chrom2 and _overlap(start1, stop1, start2, stop2) is True:
                if start1 < start2:
                    start = start1
                else:
                    start = start2
                if stop1 > stop2:
                    stop = stop1
                else:
                    stop = stop2
                reads.append(read2)
                strands.append(strand2)
                mates.append(mate2)
                d.pop(key)
                skips.append(key)
            else:
                start = start1
                stop = stop1
    try:
        start
    except NameError:  # all keys were skipped as they were merged with other reads already
        start = start1
        stop = stop1
        d[x] = [chrom1, start, stop, ref_coords, reads, strands, mates]
        skips.append(x)
    else:
        strands = list(set(strands))
        d[x] = [chrom1, start, stop, ref_coords, reads, strands, mates]
        skips.append(x)


def _overlap(start1, stop1, start2, stop2):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    for y in xrange(start2, stop2+1):
        if start1 < y < stop1:
            return True
        else:
            pass


def _get_len(infile):
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
        i, lines = _get_len(infile)
        for x in range(i):
            line = lines[x]
            field = line.rsplit()
            read1 = {'chrom': field[0], 'start': field[1], 'stop': field[2], 'strand': field[8]}
            read2 = {'chrom': field[3], 'start': field[4], 'stop': field[5], 'strand': field[9]}
            remaining = [field[10], field[11], field[12], field[13], field[14], field[15], field[6]]  # read name, TE reference coordinates, TE name, TE family
            te_coords = {'chrom': field[10], 'start': field[11], 'stop': field[12], 'strand': field[13], 'name': field[14]}
            if _overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                te_read = read1
                dna_read = read2
                mate = 1  # te read is mate 1 in paired data
            elif _overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
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


def split_bed_by_gene(bedfile, splitcol, out_prefix):
    """
    Takes input bedfile sorted by TE name and splits into
    different file for each TE name
    """
    with open(bedfile, 'r') as bedin:
        for line in bedin:
            field = line.rsplit()
            te = field[splitcol]
            try:
                prev_te
            except NameError:
                outfile = open("{pref}_{te}". format(pref=out_prefix, te=te), 'a+')
                outfile.write(line)
                prev_te = te
            else:
                if te == prev_te:
                    outfile.write(line)
                    prev_te = te
                else:
                    outfile.close()
                    outfile = open("{pref}_{te}". format(pref=out_prefix, te=te), 'a+')
                    outfile.write(line)
                    prev_te = te
        outfile.close()


def separate_reads(acc):
    """
    splits read name info into different file and adds unique IDs for insertions
    """
    with open('insertions.temp', 'r') as infile, open('insertions_unsorted.temp', 'w+') as outfile, open('id_{a}.fa'.format(a=acc), 'w+') as id_file:
        x = 0
        for line in infile:
            line = line.rsplit()
            data = line[:9]
            reads = line[9]
            mates = line[10]
            ident = acc + '_' + str(x)
            outfile.write('{data}\t{id}\n'.format(data='\t'.join(data), id=ident))
            id_file.write('>{id}\t{reads}\t{mates}\n'.format(id=ident, reads=reads, mates=mates))
            x += 1


def splitfile(inp):
    """
    splits single and double breakpoints
    """
    with open(inp, 'r') as infile, open('single_break.temp', 'w+') as single, open('double_break.temp', 'w+') as double:
        for line in infile:
            field = line.rsplit()
            count = int(field[11])
            data = field[:11]
            if count == 1:
                single.write("{data}\n".format(data='\t'.join(data)))
            elif count == 2:
                double.write("{data}\n".format(data='\t'.join(data)))
            else:
                pass


def _filter_del(inf, master, outf, ref):
    """
    Take bedfile containing all TE deletions and create
    polymorphic TE file with coordinates of TE and
    list of accessions that contain the TE
    """
    with open(inf, 'r') as infile, open(outf, 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            accessions = line[6]
            accessions = accessions.split(',')
            coords = line[:5]
            temp = [ref]
            for item in master:
                if item not in accessions:
                    temp.append(item)
                else:
                    pass
            info = '\t'.join(coords) + '\t' + ','.join(temp) + '\n'
            outfile.write(info)


def _create_names_list(inf):
    names = []
    with open(inf, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            names.append(line[0])
    return names


def inverse_del(inf, outf, ref):
    master = _create_names_list(inf)
    _filter_del(inf, master, outf, ref)


def filter_unique_break(feature):
    """
    filters where there is a breakpoint
    use with pybedtools.filter()
    """
    if int(feature[3]) == 1 and int(feature[4]) != 1:
        return True
    elif int(feature[3]) != 1 and int(feature[4]) == 1:
        return True
    else:
        return False

def break_coords(feature):
    """
    reorders split read coords with location of breakpoint
    use with pybedtools.each()
    """
    start_count = int(feature[3])
    stop_count = int(feature[4])
    start = feature[1]
    stop = feature[2]
    chrom = feature[0]
    if start_count == 1:
        feature = [chrom, start, start]
        return feature
    elif stop_count == 1:
        feature = [chrom, stop, stop]
        return feature
    else:
        raise Exception('incorrect filtering of breakpoints')


def create_deletion_coords(bedfile, saveas):  # problem where some deletion start coords are before stop
    """
    Creates set of putative deletion coordinates where discordant
    read pairs are on same chromosome, different strands, and
    are at least 3 standard deviations from the mean insert size
    and less than 20 kb from each other.
    Assumes input bedfile only contains discordant reads
    """
    with open(saveas, 'w+') as outfile:
        for line in bedfile:
            chr1 = line[0]
            start1 = int(line[1])
            stop1 = int(line[2])
            strand1 = line[8]
            chr2 = line[3]
            start2 = int(line[4])
            stop2 = int(line[5])
            read = line[6]
            strand2 = line[9]
            if chr1 == chr2 and strand1 != strand2:
                if _overlap(start1, stop1, start2, stop2) is True:
                    pass
                else:
                    if start2 >= stop1:
                        start = stop1
                        stop = start2
                    else:
                        start = stop2
                        stop = start1
                    gapsize = stop - start
                    if  gapsize < 20000:
                        outfile.write('{ch}\t{start}\t{stop}\t{read}\n'.format(ch=chr1,
                                                                               start=start,
                                                                               stop=stop,
                                                                               read=read))
                    else:
                        pass
            else:
                pass


def count_inserts(inf, outf, chrom):
    """
    Splits chromosome into 50 bins.
    Finds frequency of TE inserts in each bin.
    Returns bin start points and number of insertions in each bin as tsv file.
    """
    with open(inf, 'r') as infile, open(outf, 'w+') as outfile:
        chr_size = {1: 30427617, 2: 19698289, 3: 23459830, 4: 18585056, 5: 26975502}
        length = chr_size[chrom]
        bins = length / 50
        bins = int(bins)
        bins_dict = {}
        for x in range(50):
            bins_dict[x] = 0
        for line in infile:
            line = line.rsplit()
            ins_chr = line[0]
            ins_chr = int(ins_chr)
            ins_start = int(line[1])
            ins_end = line[2]
            if ins_chr == chrom:
                bin_no = _which_bin(bins, 50, ins_start)
                bins_dict[bin_no] += 1
            else:
                pass
        for key, value in bins_dict.items():
            outfile.write('{v}\n'.format(v=value))


def _which_bin(size, number, inp):
    """
    Takes number of bins, bin size, and a number to evaluate.
    Finds which bin the input number goes into.
    Returns bin number.
    """
    if inp <= size:
        return 1
    else:
        for x in xrange(number):
            if (x*size) <= inp <= ((x+1)*size):
                return x
            else:
                pass


def convert_split_pairbed(inp, outf):
    """
    converts split read bedfile into bedpe format
    with each read on one line
    """
    with open(inp, 'r') as infile, open(outf, 'w+') as outfile:
        i, lines = _get_len(infile)
        x = 0
        while x < i:
            coords, read, strand = _get_features(lines[x])
            x += 1
            next_coords, next_read, next_strand = _get_features(lines[x])
            if next_read == read:
                mate = read.split('_')
                outfile.write("{co}\t{nco}\t{read}\t{mt}\t{st1}\t{st2}\n".format(co='\t'.join(coords),
                                                                                 nco='\t'.join(next_coords),
                                                                                 read=mate[0],
                                                                                 mt=mate[1],
                                                                                 st1=strand,
                                                                                 st2=next_strand))
                x += 1
            else:
                pass


def _get_features(inp):
    line = inp.rsplit()
    coords = line[:3]
    read = line[3]
    strand = line[5]
    return coords, read, strand


def _get_data(inp):
    lengths = []
    for line in inp:
        length = int(line[8])
        if length > 0:
            lengths.append(length)
        else:
            pass
    return lengths


def _reject_outliers(data, m=2.):
    """
    rejects outliers more than 2
    standard deviations from the median
    """
    median = np.median(data)
    std = np.std(data)
    for item in data:
        if abs(item - median) > m * std:
            data.remove(item)
        else:
            pass


def _calc_size(data):
    mn = int(np.mean(data))
    std = int(np.std(data))
    return mn, std


def calc_mean(data):
    lengths = _get_data(data)
    _reject_outliers(lengths)
    mn, std = _calc_size(lengths)
    return mn, std


def annotate_single_breakpoint():
    """adds breakpoint coordinates to insertion"""
    with open('single_break_intersect.temp', 'r') as infile, open('insertions.temp', 'a+') as outfile:
        for line in infile:
            line = line.rsplit()
            coords = line[11:14]
            data = line[3:11]
            outfile.write('{coords}\t{data}\n'.format(coords='\t'.join(coords), data='\t'.join(data)))


def annotate_double_breakpoint():
    """adds breakpoint coordinates to insertions with two breakpoints"""
    with open('double_break_intersect.temp', 'r') as infile, open('insertions.temp', 'a+') as outfile:
        i, lines = _get_len(infile)
        x = 0
        while x < i:
            line = lines[x].rsplit()
            chrom = line[11]
            break_1 = int(line[12])
            data = line[3:11]
            x += 1
            nextline = lines[x].rsplit()
            break_2 = int(nextline[12])
            if break_1 > break_2:
                start = break_2
                stop = break_1
                outfile.write('{ch}\t{start}\t{stop}\t{data}\n'.format(ch=chrom, start=start, stop=stop, data='\t'.join(data)))
                x += 1
            elif break_1 < break_2:
                start = break_1
                stop = break_2
                outfile.write('{ch}\t{start}\t{stop}\t{data}\n'.format(ch=chrom, start=start, stop=stop, data='\t'.join(data)))
                x += 1
            else:
                raise Exception('Incorrect breakpoint information')


def annotate_deletions(inp, acc):
    """
    Calls deletions where the gap between paired reads is at
    least 75 percent the length of the TE
    """
    x = 0
    with open(inp, 'r') as infile, open('deletions_{a}.bed'.format(a=acc), 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            coords = [line[0], int(line[1]), int(line[2])]  # chr, start, stop
            te = [line[4], line[5], line[6], line[7], line[8]]  # chr, start, stop, strand, name
            overlap = int(line[11])
            gapsize = coords[2] - coords[1]
            percentage = overlap / gapsize
            if percentage >= 0.75:
                try:
                    name
                except NameError:
                    ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                    data = map(str, te)
                    outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                    x += 1
                    name = te[4]
                else:
                    if name != te[4]:
                        ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                        data = map(str, te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                        x += 1
                        name = te[4]
                    else:
                        name = te[4]
            else:
                pass


def annotate_insertions(collapse_file, insertion_file, accession_name):
    """
    Find insertion coordinates and TE orientation. Adds unique ID: <accession_name>_<number>
    having non-int chromosome names here is a problem
    assumes all read pairs are discordant
    """
    with open(collapse_file, 'r') as infile, open(insertion_file, 'w+') as outfile:
        i, lines = _get_len(infile)
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            strand = line[3]
            te_name = line[4]
            mate = line[10]
            mate = mate.split(',')
            te_reads = line[9]
            te_reads = te_reads.split(',')
            reference = [line[5], line[6], line[7], line[8]]  # reference chrom, start, stop, strand
            pair = _find_next(lines, i, x, chrom, strand, start, stop, te_name)
            if strand != reference[3]:
                if reference[3] == '+':
                    orientation = '-'
                else:
                    orientation = '+'
            else:
                orientation = reference[3]
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                pair_start = pair[0]
                pair_mates = pair[2]
                next_read_names = pair[1]
                mate = pair_mates + mate
                te_reads = next_read_names + te_reads
                outfile.write('{chr}\t{start}\t{stop}\t{orient}\t{name}\t{ref}\t{reads}\t{mates}\n'.format(chr=chrom,
                                                                                                           start=stop,
                                                                                                           stop=pair_start,
                                                                                                           orient=orientation,
                                                                                                           name=te_name,
                                                                                                           ref='\t'.join(reference),
                                                                                                           reads='|'.join(te_reads),
                                                                                                           mates='|'.join(mate)))


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
def _find_next(lines, i, x, chrom, strand, start, stop, te_name):
    """
    Find next read linked to same TE. Looks in 100 bp window.
    """
    while True:
        line = lines[x+1]
        line = line.rsplit()
        next_chrom = line[0]
        next_start = int(line[1])
        next_stop = int(line[2])
        next_strand = line[3]
        next_te_name = line[4]
        next_mate = line[10]
        next_mate = next_mate.split(',')
        next_te_reads = line[9]
        next_te_reads = next_te_reads.split(',')
        if strand != next_strand and te_name == next_te_name and chrom == next_chrom and stop <= next_start and (stop + 100) > next_start:
            return next_start, next_te_reads, next_mate
        elif stop + 100 < next_start:
            return False
        else:
            x += 1
            if x >= i:
                return False
