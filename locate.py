from __future__ import division
import os
from sys import argv
import numpy as np
import pysam


def _overlap(start1, stop1, start2, stop2, d=0):
    """
    Returns True if sets of coordinates overlap.
    Assumes coordinates are on same chromosome.
    10 bp window (seems to work better)
    """
    d = int(d)
    for y in xrange(start2-d, stop2+d):
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
            r1 = True if _overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop']), 10) is True else False
            r2 = True if _overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop']), 10) is True else False
            if r1 is True and r2 is False:
                dna_read = read2
                te_read = [read1['chrom'], read1['start'], read1['stop']]
                mate = 2
            elif r1 is False and r2 is True:
                dna_read = read1
                te_read = [read2['chrom'], read2['start'], read2['stop']]
                mate = 1
            else:
                pass
            # bedpe format demands chr-start-stop-chr-start-stop-strand1-strand2
            write_string = '{chr1}\t{start1}\t{stop1}\t{chr2}\t{start2}\t{stop2}\t{strand1}\t{strand2}\t{rd}\t{te}\t{sd}\t{te_read}\n'.format(
                chr1=dna_read['chrom'],
                start1=dna_read['start'],
                stop1=dna_read['stop'],
                chr2=te_coords['chrom'],
                start2=te_coords['start'],
                stop2=te_coords['stop'],
                strand1=dna_read['strand'],
                strand2=te_coords['strand'],
                rd=field[6],
                te=field[15],
                sd=sd,
                te_read='\t'.join(te_read))
            if sd == 'disc':
                if (mate == 1 and dna_read['strand'] == '+') or (mate == 2 and dna_read['strand'] == '-'):
                    disc_forward.write(write_string)
                elif (mate == 1 and dna_read['strand'] == '-') or (mate == 2 and dna_read['strand'] == '+'):
                    disc_reverse.write(write_string)
            else:
                split_out.write(write_string)


def _condense_coords(starts, stops):
    """
    finds if lists of coordinate sets all overlap one another
    Returns widest set of merged coordinates (to encompass all overlapping TEs)
    Returns None if there is no overlap
    """
    start = starts[0]
    stop = stops[0]
    for x in range(1, len(starts)):
        if _overlap(start, stop, starts[x], stops[x], 0) is True:
            if start > starts[x]:
                start = starts[x]
            else:
                continue
            if stops[x] < stops[x]:
                stop = stops[x]
            else:
                continue
        else:
            break
    else:
        return start, stop


def process_merged(infile, outfile, sd):
    """
    take merged coordinates and filter out those where multiple non-nested TEs insert into same locus
    resolve cases where insertion site is within another TE
    """
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf:
        for line in inf:
            line = line.rsplit()
            te_chroms = line[3].split(',')
            te_names = line[7].split(',')
            r2_chrom = line[9].split(',')
            r2_start = line[10].split(',')
            r2Starts = [int(x) for x in r2_start]
            r2_stop = line[11].split(',')
            r2Stops = [int(x) for x in r2_stop]
            if len(te_chroms) > 1 or len(r2_chrom) > 1:
                pass
            else:
                r2_coords = _condense_coords(r2Starts, r2Stops)
                starts = line[4].split(',')
                starts = [int(x) for x in starts]
                stops = line[5].split(',')
                stops = [int(x) for x in stops]
                if r2_coords is None or sd == 'disc' or (abs(max(stops) - min(starts) <= 200)):
                    if len(te_names) > 1:
                        coords = _condense_coords(starts, stops)
                        if coords is not None:
                            outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{cnt}\t{sd}\n'.format(ch=line[0],
                                                                                                                sta=line[1],
                                                                                                                stp=line[2],
                                                                                                                tec=line[3],
                                                                                                                tesa=coords[0],
                                                                                                                tesp=coords[1],
                                                                                                                rds=line[6],
                                                                                                                nm=line[7],
                                                                                                                cnt=line[8],
                                                                                                                sd=sd))
                        else:
                            pass
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
                                                                                                                sd=sd))
                else:
                    pass


def process_merged_disc(infile, outfile):
    """
    takes merged coordinates and finds where there are discordant reads in both direction
    collects read count information and writes to file
    """
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf:
        for line in inf:
            line = line.rsplit()
            te_chroms = line[3].split(',')
            te_names = line[7].split(',')
            read_types = line[9].split(',')
            starts = line[4].split(',')
            stops = line[5].split(',')
            if len(te_chroms) > 1:
                pass
            elif len(te_names) > 1:
                starts = [int(x) for x in starts]
                stops = [int(x) for x in stops]
                coords = _condense_coords(starts, stops)
                if coords is not None and len(read_types) == 2:
                    outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{count}\n'.format(
                        ch=line[0],
                        sta=line[1],
                        stp=line[2],
                        tec=line[3],
                        tesa=coords[0],
                        tesp=coords[1],
                        rds=line[6],
                        nm=line[7],
                        count=line[8]))
                else:
                    pass
            elif len(read_types) == 2:
                start = starts[0]
                stop = stops[0]
                outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{count}\n'.format(
                        ch=line[0],
                        sta=line[1],
                        stp=line[2],
                        tec=line[3],
                        tesa=start,
                        tesp=stop,
                        rds=line[6],
                        nm=line[7],
                        count=line[8]))
            else:
                pass


def separate_reads(infile, outfile, reads_file):
    """
    splits read name info into different file and adds unique IDs for insertions
    """
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf, open(reads_file, 'w+') as fasta_file:
        x = 0
        for line in inf:
            line = line.rsplit()
            data = line[:6]
            name = line[7]
            reads = line[6]
            outf.write('{data}\t{te}\t{x}\n'.format(data='\t'.join(data), te=name, x=x))
            fasta_file.write('>{x}\t{reads}\n'.format(x=x, reads=reads))
            x += 1


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
            if chr1 == chr2:
                if _overlap(start1, stop1, start2, stop2, 10) is True:
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


def calc_cov(bam_name, start, stop):
    """
    calculates average coverage
    """
    bam = pysam.AlignmentFile(bam_name)
    # get chromosome names
    nms = []
    for i in bam.header['SQ']:
        nms.append(i['SN'])
    x = 0
    l = 0
    for read in bam.pileup(nms[0], start, stop):
        x += read.n
        l += 1
    if l > 0 and x > 0:
        bam.close()
        return int(x/l)
    else:
        bam.close()
        return 0


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


def annotate_deletions(inp, acc, num_reads, bam, mn):
    """
    Calls deletions where the gap between paired reads is at
    least 20 percent the length of the TE
    and there are either:
       1 split/disc read spanning the TE and
       coverage at TE is 1/10 the coverage in surrounding area, or
       num_reads split reads spanning the TE
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
                tes[name] = [cov, 0, 0]  # coverage, split, disc
            else:
                pass
            if (gapsize <= 0) or (name in written_tes) or ((length-mn) > gapsize):
                pass
            else:
                percentage = overlap / gapsize
                if (read_type == 'split' and percentage >= 0.7) or (read_type == 'disc' and percentage >= 0.2):
                    if read_type == 'split':
                        tes[name][1] += 1
                    elif read_type == 'disc':
                        tes[name][2] += 0.5
                    else:
                        raise Exception('Incorrect read type information')
                    split = tes[name][1]
                    disc = tes[name][2]
                    total_reads = split + disc
                    if tes[name][1] > 0 and tes[name][2] > 0:
                        both_types = True
                    else:
                        both_types = False
                    if tes[name][0] <= 0.1 or (total_reads >= num_reads) or (length <= 1000 and total_reads >= (num_reads/2)): #or (both_types == True and total_reads >= (num_reads/2)):
                        ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                        data = (str(x) for x in te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                        x += 1
                        written_tes.append(name)
                    else:
                        pass
                else:
                    pass
    test_head.close()
    allreads.close()


def append_origin(feature, word):
    """
    use with pybedtools.each()
    append 'word' as final column in file
    """
    feature.append(word)
    return feature


def condense_names(feature):
    """
    use in pybedtools.each()
    """
    feature = feature[:-1]
    names = set(feature[-1].split(','))
    names = ','.join(names)
    feature[-1] = names
    return feature

def reorder_intersections(feature, num_disc, num_split):
    """
    use with pybedtools.each()
    """
    chrom = feature[0]
    start = feature[1]
    stop = feature[2]
    techrom = feature[13]
    testart = feature[14]
    testop = feature[15]
    reads = set(feature[6].split(',') + feature[16].split(','))
    names = set(feature[7].split(',') + feature[17].split(','))
    disc_reads = feature[8]
    split_reads = feature[-2]
    if disc_reads >= num_disc and split_reads >= num_split:
        feature = [chrom, start, stop, techrom, testart, testop, ','.join(reads), ','.join(names)]
        return feature
    else:
        pass
