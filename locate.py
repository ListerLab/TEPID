from __future__ import division
import os
from sys import argv
import numpy as np
import pysam


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

# def append_break(feature):
#     """
#     use with pybedtools.each()
#     """
#     if int(feature[-1]) >= 1:
#         breakpoint = 'True'
#     else:
#         breakpoint = 'False'
#     feature = feature[:11]
#     feature.append(breakpoint)
#     return feature


def _overlap(start1, stop1, start2, stop2):
    """
    Returns True if sets of coordinates overlap.
    Assumes coordinates are on same chromosome.
    10 bp window (seems to work better)
    """
    for y in xrange(start2-10, stop2+10):
        if start1 <= y <= stop1:
            return True
        else:
            pass


def _get_len(infile):
    """returns number of lines in file and all lines as part of list"""
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    try:
        return i, lines
    except:
        return 0, 0


def reorder(insert_file, split_outf, disc_forw, disc_rev):
    """
    Reorder columns so that TE read is in second position.
    """
    with open(insert_file, 'r') as infile, open(split_outf, 'w+') as split_out, open(disc_forw, 'w+') as disc_forward, open(disc_rev, 'w+') as disc_reverse:
        for line in infile:
            field = line.rsplit()
            read1 = {'chrom': field[0], 'start': field[1], 'stop': field[2], 'strand': field[8]}
            read2 = {'chrom': field[3], 'start': field[4], 'stop': field[5], 'strand': field[9]}
            sd = field[10]
            te_coords = {'chrom': field[11], 'start': field[12], 'stop': field[13], 'strand': field[14], 'name': field[15]}
            if _overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                dna_read = read2
                mate = 2  # DNA read is mate 1 in paired data
            elif _overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                dna_read = read1
                mate = 1
            else:
                raise Exception('check coords')
            # bedpe format demands chr-start-stop-chr-start-stop-strand1-strand2
            write_string = '{chr1}\t{start1}\t{stop1}\t{chr2}\t{start2}\t{stop2}\t{strand1}\t{strand2}\t{rd}\t{te}\t{sd}\n'.format(chr1=dna_read['chrom'],
                                                                                                                          start1=dna_read['start'],
                                                                                                                          stop1=dna_read['stop'],
                                                                                                                          chr2=te_coords['chrom'],
                                                                                                                          start2=te_coords['start'],
                                                                                                                          stop2=te_coords['stop'],
                                                                                                                          strand1=dna_read['strand'],
                                                                                                                          strand2=te_coords['strand'],
                                                                                                                          rd=field[6],
                                                                                                                          te=field[15],
                                                                                                                          sd=sd)
            if sd == 'disc':
                if (mate == 1 and dna_read['strand'] == '+') or (mate == 2 and dna_read['strand'] == '-'):
                    disc_forward.write(write_string)
                elif (mate == 1 and dna_read['strand'] == '-') or (mate == 2 and dna_read['strand'] == '+'):
                    disc_reverse.write(write_string)
            else:
                split_out.write(write_string)


def process_merged(infile, outfile, sd):
    """
    take merged coordinates and filter out those where multiple non-nested TEs insert into same locus
    """
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf:
        for line in inf:
            line = line.rsplit()
            te_chroms = line[3].split(',')
            te_names = line[7].split(',')
            if len(te_chroms) > 1:
                pass
            elif len(te_names) > 1:
                starts = line[4].split(',')
                starts = [int(x) for x in starts]
                stops = line[5].split(',')
                stops = [int(x) for x in stops]
                start = starts[0]
                stop = stops[0]
                for x in range(1, len(starts)):
                    if _overlap(start, stop, starts[x], stops[x]) is True:
                        if start < starts[x]:
                            start = starts[x]
                        else:
                            continue
                        if stops[x] > stops[x]:
                            stop = stops[x]
                        else:
                            continue
                    else:
                        break
                else:
                    outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{cnt}\t{sd}\n'.format(ch=line[0],
                                                                                                        sta=line[1],
                                                                                                        stp=line[2],
                                                                                                        tec=line[3],
                                                                                                        tesa=start,
                                                                                                        tesp=stop,
                                                                                                        rds=line[6],
                                                                                                        nm=line[7],
                                                                                                        cnt=line[8],
                                                                                                        sd=sd
                                                                                                        ))
            else:
                start = line[4].split(',')
                start = start[0]
                stop = line[5].split(',')
                stop = stop[0]
                outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{cnt}\t{sd}\n'.format(ch=line[0],
                                                                                                        sta=line[1],
                                                                                                        stp=line[2],
                                                                                                        tec=line[3],
                                                                                                        tesa=start,
                                                                                                        tesp=stop,
                                                                                                        rds=line[6],
                                                                                                        nm=line[7],
                                                                                                        cnt=line[8],
                                                                                                        sd=sd
                                                                                                        ))



def separate_reads(acc):
    """
    splits read name info into different file and adds unique IDs for insertions
    """
    with open('insertions.temp', 'r') as infile, open('insertions_unsorted.temp', 'w+') as outfile, open('id_{a}.fa'.format(a=acc), 'w+') as id_file:
        x = 0
        for line in infile:
            line = line.rsplit()
            data = line[:8]
            reads = line[8]
            mates = line[9]
            ident = acc + '_' + str(x)
            outfile.write('{data}\t{id}\n'.format(data='\t'.join(data), id=ident))
            id_file.write('>{id}\t{reads}\t{mates}\n'.format(id=ident, reads=reads, mates=mates))
            x += 1


# def filter_unique_break(feature):
#     """
#     filters where there is a breakpoint
#     use with pybedtools.filter()
#     """
#     if int(feature[3]) == 1 and int(feature[4]) != 1:
#         return True
#     elif int(feature[3]) != 1 and int(feature[4]) == 1:
#         return True
#     else:
#         return False


def filter_discordant(bam, dist, new_filename):
    """
    filters discordant reads in bamfile and writes to new bam
    """
    bamfile = pysam.AlignmentFile(bam, 'rb')
    header = bamfile.header.copy()
    new_bam = pysam.Samfile(new_filename, 'wb', header=header)
    for i in bamfile:
        if abs(i.tlen) > dist or i.reference_id != i.next_reference_id:
            new_bam.write(i)
        else:
            pass
    new_bam.close()
    bamfile.close()


# def break_coords(feature):
#     """
#     reorders split read coords with location of breakpoint
#     use with pybedtools.each()
#     """
#     start_count = int(feature[3])
#     stop_count = int(feature[4])
#     start = feature[1]
#     stop = feature[2]
#     chrom = feature[0]
#     if start_count == 1:
#         feature = [chrom, start, start]
#         return feature
#     elif stop_count == 1:
#         feature = [chrom, stop, stop]
#         return feature
#     else:
#         raise Exception('incorrect filtering of breakpoints')


def create_deletion_coords(bedfile, saveas):
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
            read_type = line[10]
            # could add the breakpoint info
            if chr1 == chr2:
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
                        outfile.write('{ch}\t{start}\t{stop}\t{read}\t{rt}\n'.format(ch=chr1,
                                                                                     start=start,
                                                                                     stop=stop,
                                                                                     read=read,
                                                                                     rt=read_type))
                    else:
                        pass
            else:
                pass


def convert_split_pairbed(inp, outf):
    """
    converts split read bedfile into bedpe format
    with each read on one line
    read names need to be order
    """
    with open(inp, 'r') as infile, open(outf, 'w+') as outfile:
        i, lines = _get_len(infile)
        x = 0
        while x < i:
            coords, read, strand = _get_features(lines[x])
            x += 1
            next_coords, next_read, next_strand = _get_features(lines[x])
            if next_read == read:
                mate = read[-1]
                rd = read[:-2]
                outfile.write("{co}\t{nco}\t{read}\t{mt}\t{st1}\t{st2}\n".format(co='\t'.join(coords),
                                                                                 nco='\t'.join(next_coords),
                                                                                 read=rd,
                                                                                 mt=mate,
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


def get_coverages(chrom, start, stop, bam, chrom_sizes):
    """
    find average coverage in given region
    compared to +/- 2kb surrounding region
    """
    te = 0
    l = 0
    ustream = 0
    ul = 0
    dstream = 0
    dl = 0

    if (start - 2000) > 0:
        ustart = (start - 2000)
    else:
        ustart = 0

    if (stop + 2000) < chrom_sizes[chrom]:
        dstop = stop + 2000
    else:
        dstop = chrom_sizes[chrom]

    for read in bam.pileup(chrom, start, stop):
        te += read.n
        l += 1
    for read in bam.pileup(chrom, ustart, start):
        ustream += read.n
        ul += 1
    for read in bam.pileup(chrom, stop, dstop):
        dstream += read.n
        dl += 1
    if (ustream + dstream) > 0:
        surround = (ustream + dstream) / (ul + dl)
    else:
        ratio = 0
    if te > 0:
        tot_te = te / l
        ratio =  tot_te / surround
    else:
        ratio = 0
    return ratio


def annotate_deletions(inp, acc, num_split, bam, mn):
    """
    Calls deletions where the gap between paired reads is at
    least 20 percent the length of the TE
    and there are either:
       one discordant read pair spanning TE, or
       1 split read spanning the TE and
       coverage at TE is 1/10 the coverage in surrounding area, or
       num_split split reads spanning the TE
    """
    x = 0
    tes = {}
    written_tes = []

    # check if sorted
    test_head = pysam.AlignmentFile(bam, 'rb')
    chrom_sizes = {}
    for i in test_head.header['SQ']:
        chrom_sizes[i['SN']] = int(i['LN'])
    if test_head.header['HD']['SO'] == 'coordinate':
        pass
    else:
        print 'Sorting bam file'
        pysam.sort(bam, 'sorted.temp')
        os.remove(bam)
        os.rename('sorted.temp.bam', bam)
    
    # check if indexed
    if '{}.bai'.format(bam) in os.listdir('.'):
        print '  Using index {}.bai'.format(bam)
        allreads = pysam.AlignmentFile(bam, 'rb')
    else:
        print '  Indexing bam file'
        pysam.index(bam)
        allreads = pysam.AlignmentFile(bam, 'rb')

    with open(inp, 'r') as infile, open('deletions_{a}.bed'.format(a=acc), 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            coords = [line[0], int(line[1]), int(line[2])]  # chr, start, stop
            te = [line[5], line[6], line[7], line[8], line[9]]  # chr, start, stop, strand, name
            name = te[4]
            length = int(te[2]) - int(te[1])
            overlap = int(line[12])
            gapsize = coords[2] - coords[1]
            read_type = line[4]
            if name not in tes.keys():
                cov = get_coverages(coords[0], coords[1], coords[2], allreads, chrom_sizes)
                tes[name] = [cov, 0]
            else:
                pass
            if (gapsize <= 0) or (name in written_tes) or ((length-mn) > gapsize):
                pass
            else:
                percentage = overlap / gapsize
                if percentage >= 0.2:  # 0.2 best so far
                    tes[name][1] += 1
                    if (tes[name][0] <= 0.1) or (tes[name][1] >= num_split) or (length <= 1000 and tes[name][1] >= (num_split/2)) or (read_type == 'disc'):
                        ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                        data = (str(x) for x in te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                        x += 1
                        written_tes.append(name)
                    else:
                        pass
                else:
                    pass


def append_origin(feature, word):
    """
    use with pybedtools.each()
    append 'word' as final column in file
    """
    feature.append(word)
    return feature


def annotate_insertions(collapse_file, insertion_file, id_file, num_reads):
    """
    Find insertion coordinates and TE orientation.
    assumes all read pairs are discordant
    """
    with open(collapse_file, 'r') as infile, open(insertion_file, 'w+') as outfile, open(id_file, 'w+') as ids:
        i, lines = _get_len(infile)
        if i > 0:
            for x in range(i):
                line = lines[x]
                line = line.rsplit()
                chrom = line[0]
                start = int(line[1])
                stop = int(line[2])
                te_chrom = line[3]
                te_start = line[4]
                te_stop = line[5]
                te_reads = line[6]
                te_name = line[7].split(',')
                count = int(line[8])
                sd = line[9]
                pair = _find_next(lines, i, x, chrom, start, stop, te_name, sd)
                if pair is False:
                    if len(sd) > num_reads:
                        outfile.write('{chr}\t{start}\t{stop}\t{name}\t{te_chr}\t{te_start}\t{te_stop}\n'.format(chr=chrom,
                                                                                                                 start=start,
                                                                                                                 stop=stop,
                                                                                                                 name=','.join(te_name),
                                                                                                                 te_chr=te_chrom,
                                                                                                                 te_start=te_start,
                                                                                                                 te_stop=te_stop))
                else:
                    pair_start = int(pair[0])
                    if pair_start < stop:
                        pair_start = stop
                    next_read_names = pair[1]
                    te_reads = next_read_names + te_reads
                    sd = sd + pair[2]
                    outfile.write('{chr}\t{start}\t{stop}\t{name}\t{te_chr}\t{te_start}\t{te_stop}\n'.format(chr=chrom,
                                                                                                             start=stop,
                                                                                                             stop=pair_start,
                                                                                                             name=','.join(te_name),
                                                                                                             te_chr=te_chrom,
                                                                                                             te_start=te_start,
                                                                                                             te_stop=te_stop))
        else:
            pass

 
def _find_next(lines, i, x, chrom, start, stop, te_name, sd):
    """
    Find next read linked to same TE. Looks in 10 bp window.
    """
    if x+1 < i:
        while True:
            line = lines[x+1]
            line = line.rsplit()
            next_chrom = line[0]
            next_start = int(line[1])
            next_stop = int(line[2])
            next_te_chrom = line[3]
            next_te_start = line[4]
            next_te_stop = line[5]
            next_te_reads = line[6]
            next_te_name = line[7].split(',')
            next_count = int(line[8])
            next_sd = line[9]
            correct_te = False
            for n in next_te_name:
                if n in te_name:
                    correct_te = True
                else:
                    pass
            if correct_te is True and chrom == next_chrom and _overlap(start, stop, next_start, next_stop) is True:
                return next_start, next_te_reads, next_sd
            elif stop + 100 < next_start:
                return False
            else:
                x += 1
                if x >= i:
                    return False
    else:
        return False
