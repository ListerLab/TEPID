#!/usr/bin/env python

from __future__ import division
import os
import numpy as np
import pysam
import heapq
import pybedtools
from glob import glob
from time import ctime


def readNames(names):
    n = []
    with open(names, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            if len(line) > 0:
                n.append(line[0])
            else:
                pass
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


def flatten_list(l):
    s = []
    for i in l:
        if i not in s:
            s.append(i)
    return s


def check_no_breaks(coords):
    """[[start,stop]...[start,stop]]"""
    copy = sorted(list(coords))
    for x in range(len(copy)-1):
        if _overlap(copy[x][0], copy[x][1], copy[x+1][0], copy[x+1][1], 0) is False:
            return False
        else:
            pass
    return True


def find_reads_conc(coords, bam):
    """
    If there are concordant reads spanning input coordinates return true, else return False
    """
    reads = bam.fetch(coords[0], coords[1], coords[2])
    intervals = []
    for read in reads:
        if read.pos is not None and read.aend is not None:
            intervals.append([read.pos, read.aend])
    intervals = flatten_list(intervals)
    return check_no_breaks(intervals)


def extract_reads(bam, name_indexed, names, acc):
    header = bam.header.copy()
    out_name = 'extracted_reads_{}.bam'.format(acc)
    out_bam = pysam.Samfile(out_name, 'wb', header=header)
    for name in names:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out_bam.write(x)
    out_bam.close()
    return out_name


def check_te_overlaps(te, bamfile, te_list):
    intersections = pybedtools.BedTool(bamfile).bam_to_bed().intersect(te, wb=True, nonamecheck=True)
    reads = []
    for r in intersections:
        if r[-3] in te_list:
            reads.append(r[3])
        else:
            pass
    return reads


def check_te_overlaps_dels(te, bamfile, te_list):
    pybedtools.BedTool(bamfile).bam_to_bed().saveas('split.temp')
    convert_split_pairbed('split.temp', 'split_bedpe.temp')
    split = pybedtools.BedTool('split_bedpe.temp').each(append_origin, word='split').saveas()
    create_deletion_coords(split, 'second_pass_del_coords.temp')
    dels = pybedtools.BedTool('second_pass_del_coords.temp').sort()
    intersections = dels.intersect(te, wb=True, nonamecheck=True, sorted=True)
    os.remove('second_pass_del_coords.temp')
    reads = []
    for r in intersections:
        if r[-3] in te_list:
            reads.append(r[3])
        else:
            pass
    return reads


def get_last_id(acc, indel):
    with open("{i}_reads_{a}.txt".format(i=indel, a=acc), 'r') as f:
        for line in f:
            pass
    return int(line.rsplit()[0][1:])


def write_te(te_file, read_file, data, read_names, iterator):
    # chr start stop chr start stop te accessions
    coords = [str(x) for x in data[:-2]]
    te_file.write("\t".join(coords)+"\t"+data[-2]+"\t"+str(iterator)+"\n")
    read_file.write(">"+str(iterator)+"\t"+",".join(read_names)+"\n")


def ambiguous(coords, concordant):
    """
    check if there is not enough evidence to give confident call
    return True if call is ambiguous
    """
    # calculate coverage over region
    read_count = 0
    chrom = coords[0]
    start = coords[1]-200
    stop = coords[2]+200
    length = stop-start
    for read in bam.pileup(chrom, start, stop):
        read_count += read.n
    av = read_count / length
    # if lower than threshold, add to ambiguous calls
    if av < 8:
        return True
    else:
        return False


def write_ambiguous(fname, data):
    coords = [str(x) for x in data[:-2]]
    te_file.write("\t".join(coords)+"\t"+data[-2]+"\n")


def process_missed(data, indel, concordant, split_alignments, name_indexed, acc, te, refine_read_count):
    read_file_name = "second_pass_reads_{t}_{a}.txt".format(t=indel, a=acc)
    te_file_name = "second_pass_{t}_{a}.bed".format(t=indel, a=acc)
    ambiguous_file_name = 'ambiguous_calls_{a}.bed'.format(a=acc)
    with open(read_file_name, 'w+') as read_file, open(te_file_name, 'w+') as te_file, open(ambiguous_file_name, 'w+') as no_call:
        for i in data:
            coords = (i[0][0], int(i[0][1]), int(i[0][2]))
            te_list = i[0][-2].split(",")
            if acc in i[1]:
                if find_reads_conc(coords, concordant) is False:  # returns False if there is a break in coverage within region
                    split_names = find_reads(coords, split_alignments)
                    if split_names is not False:
                        extracted = extract_reads(split_alignments, name_indexed, split_names, acc)
                        if indel == 'insertion':
                            read_names = check_te_overlaps(te, extracted, te_list)
                        elif indel == 'deletion':
                            read_names = check_te_overlaps_dels(te, extracted, te_list)
                        else:
                            raise Exception()
                        if len(read_names) >= refine_read_count:  # enough evidence to call indel
                            try:
                                iterator
                            except NameError:
                                iterator = get_last_id(acc, indel)
                            else:
                                iterator += 1
                            write_te(te_file, read_file, i[0], read_names, iterator)
                        # not enough split read to call indel, check for ambiguous call
                        elif ambiguous(coords, concordant) is True:
                            write_ambiguous(no_call, i[0])
                        else:
                            pass
                    # no split reads
                    elif ambiguous(coords, concordant) is True:
                        write_ambiguous(no_call, i[0])
                    else:
                        pass
                else:
                    pass
            else:
                pass


def refine(options):
    """
    Refine TE insertion and deletion calls within a group of related samples
    Use indel calls from other samples in the group, inspect areas of the genome in samples where
    indel was not called, and look for evidence of the same indel with much lower
    read count threshold
    """
    te = pybedtools.BedTool(options.te).sort()
    names = readNames(options.all_samples)
    if options.insertions is not False:
        insertions = getOtherLines(names, options.insertions)
    if options.deletions is not False:
        deletions = getOtherLines(names, options.deletions)  # format ([data], [inverse_accessions])
    for acc in os.listdir('.'):
        if os.path.isdir(acc) is True:
            os.chdir(acc)
            if os.path.isfile('deletions_{}.bed'.format(acc)) is True:
                print "Processing "+acc
                conc = acc+".bam"
                split = acc+".split.bam"
                check_bam(conc, options.proc)
                check_bam(split, options.proc, make_new_index=True)
                cov = calc_cov(conc, 100000, 120000)
                concordant = pysam.AlignmentFile(conc, 'rb')
                split_alignments = pysam.AlignmentFile(split, 'rb')
                name_indexed = pysam.IndexedReads(split_alignments)
                name_indexed.build()
                if options.deletions is not False:
                    print "  deletions"
                    process_missed(deletions, "deletion", concordant, split_alignments, name_indexed, acc, te, cov/5)
                else:
                    pass
                if options.insertions is not False:
                    print "  insertions"
                    process_missed(insertions, "insertion", concordant, split_alignments, name_indexed, acc, te, cov/10)
                else:
                    pass
                os.chdir('..')
            else:
                os.chdir('..')
        else:
            pass


def _overlap(start1, stop1, start2, stop2, d=0):
    """
    Returns True if sets of coordinates overlap.  Assumes coordinates are on same chromosome.
    set window size with d
    """
    d = int(d)
    for y in xrange(start2-d, stop2+d):
        if start1 <= y <= stop1:
            return True
        else:
            pass
    return False


def _get_len(infile):
    """
    returns number of lines in file and all lines as part of list
    """
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
            r1 = _overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop']), 10)
            r2 = _overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop']), 10)
            if r1 is True and r2 is False:  # find which read overlaps TE annotation
                dna_read = read2
                te_read = [read1['chrom'], read1['start'], read1['stop']]
                mate = 2
            elif r1 is False and r2 is True:
                dna_read = read1
                te_read = [read2['chrom'], read2['start'], read2['stop']]
                mate = 1
            elif r1 is True and r2 is True:
                continue  # don't write to file...need to check if this works (go to next item in for loop)
            else:
                raise Exception("Error in read cluster organization")
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


def _condense_coords(coords, d):
    """
    Merge sets of coordinates and return list of clusters with read counts for each cluster
    Input is list of coordinate tuples sorted by chrom, start
    """
    d = int(d)
    clusters = []  # populate list with individual cluster coordinates
    c = 1  # count number of reads in each cluster

    #initialize values
    chr1 = coords[0][0]
    st1 = coords[0][1]
    sto1 = coords[0][2]
    if len(coords[0]) == 4:
        te1 = [coords[0][3]]
        te = True
    else:
        te = False

    for x in range(1, len(coords)):
        # unpack values
        chr2 = coords[x][0]
        st2 = coords[x][1]
        sto2 = coords[x][2]
        if te is True:
            te2 = coords[x][3]
        else:
            pass
        if chr1 == chr2 and _overlap(st1, sto1, st2, sto2, d) is True:
            c += 1
            if te is True and te2 not in te1:
                te1.append(te2)
            else:
                pass
            if st1 > st2:
                st1 = st2
            else:
                continue
            if sto1 < sto2:
                sto1 = sto2
            else:
                continue
        else:
            # save cluster coords and reset counter
            if te is True:
                clusters.append([chr1, st1, sto1, c, te1])
                te1 = [te2]
            else:
                clusters.append([chr1, st1, sto1, c])
            c = 1
            st1 = st2
            sto1 = sto2
            chr1 = chr2
    else:
        # write final cluster
        if te is True:
            clusters.append([chr1, st1, sto1, c, te1])
        else:
            clusters.append([chr1, st1, sto1, c])
        return clusters


def process_merged(infile, outfile, sd):
    """
    take merged coordinates and filter out those where multiple non-nested TEs insert into same locus
    resolve cases where insertion site is within another TE
    """
    disc = ['disc_forward', 'disc_reverse']
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf:
        for line in inf:
            line = line.rsplit()
            te_chroms = line[3].split(',')
            te_names = line[7].split(',')
            r2_chrom = line[9].split(',')  # TE read, at end of file. Should be two clusers for split reads if TE large enough
            r2_start = line[10].split(',')
            r2Starts = [int(x) for x in r2_start]
            r2_stop = line[11].split(',')
            r2Stops = [int(x) for x in r2_stop]

            # create list of tuples to record read positions
            r2 = []
            for x in xrange(len(r2_chrom)):
                r2.append((r2_chrom[x], r2Starts[x], r2Stops[x], te_names[x]))
            # now sort reads by chr, start
            r2 = sorted(r2, key=lambda x: (x[0], x[1]))
            r2_coords = _condense_coords(r2, 10)

            starts = line[4].split(',')
            starts = [int(x) for x in starts]
            stops = line[5].split(',')
            stops = [int(x) for x in stops]
            te_reads = []
            for x in xrange(len(te_chroms)):
                te_reads.append((te_chroms[x], starts[x], stops[x], te_names[x]))
            te_reads = sorted(te_reads, key=lambda x: (x[0], x[1]))

            # filter where there are at least two clusters of reads in the TE if the TE length is more than 200 bp and read is split
            if (len(r2_coords) >= 2) or sd in disc or (abs(max(stops) - min(starts) <= 200)):
                coords = _condense_coords(te_reads, 1000)
                if len(coords) > 1:
                    crd, max_reads, total_smaller_clusters, dubs = get_main_cluster(coords)
                    te_name = ','.join(crd[3])
                else:
                    crd = coords[0][0:3]
                    max_reads = coords[0][3]
                    te_name = ','.join(coords[0][4])
                    total_smaller_clusters = 0
                    dubs = 0
                if total_smaller_clusters >= (max_reads - 2) or dubs >= 2:
                    pass
                else:
                    outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{cnt}\t{sd}\n'.format(ch=line[0],
                                                                                                        sta=line[1],
                                                                                                        stp=line[2],
                                                                                                        tec=crd[0],
                                                                                                        tesa=crd[1],
                                                                                                        tesp=crd[2],
                                                                                                        rds=line[6],
                                                                                                        nm=te_name,
                                                                                                        cnt=line[8],
                                                                                                        sd=sd))
            else:
                pass


def get_main_cluster(clusters):
    """
    takes list of read cluster coordinates and return biggest cluster,
    whether there are any of equal size,
    and whether there are smaller clusters above a certain threshold
    """
    max_read, second_max = heapq.nlargest(2, (i[3] for i in clusters))
    dubs = 0
    for i in clusters:
        if i[3] != max_read:
            pass
        else:
            dubs += 1
            chrom = i[0]
            start = i[1]
            stop = i[2]
            te = i[4]
            coords = [chrom, start, stop, te]
    return coords, max_read, second_max, dubs


def process_merged_disc(infile, outfile, num_reads, max_dist, rd_len):
    """
    takes merged coordinates and finds where there are discordant reads in both direction, linked to same TE
    collects read count information and writes to file
    filters out clusters of reads that span a distance greater than 2*(max_dist+rd_len)
    """
    with open(infile, 'r') as inf, open(outfile, 'w+') as outf:
        for line in inf:
            line = line.rsplit()
            span = int(line[2]) - int(line[1])
            max_span = (2 * (max_dist+rd_len))
            te_chroms = line[3].split(',')
            te_names = line[7].split(',')
            read_count = int(line[8])
            read_types = line[9].split(',')
            starts = line[4].split(',')
            stops = line[5].split(',')
            starts = [int(x) for x in starts]
            stops = [int(x) for x in stops]
            te_reads = []
            for x in xrange(len(te_chroms)):
                te_reads.append((te_chroms[x], starts[x], stops[x]))
            te_reads = sorted(te_reads, key=lambda x: (x[0], x[1]))
            coords = _condense_coords(te_reads, 1000)
            if len(coords) > 1:  # multiple non-overlapping TEs
                pass
            else:
                crd = coords[0][0:3]
                if read_count >= num_reads and span < max_span:
                    outf.write('{ch}\t{sta}\t{stp}\t{tec}\t{tesa}\t{tesp}\t{rds}\t{nm}\t{count}\n'.format(
                        ch=line[0],
                        sta=line[1],
                        stp=line[2],
                        tec=crd[0],
                        tesa=crd[1],
                        tesp=crd[2],
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
    are at least 4 standard deviations from the mean insert size
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
    must be sorted by read name
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


def calc_mean(bam_name, p):
    """
    calculates mean and standard deviation insert size
    """
    check_bam(bam_name, p)
    bam = pysam.AlignmentFile(bam_name)
    lengths = []
    rd = []
    x = 0
    for i in bam:
        if 10000 > i.tlen > 0 and i.tid == i.mrnm:
            lengths.append(i.tlen)
        else:
            pass
        x += 1
        rd.append(i.qlen)
        if x > 20000:
            break
    bam.close()
    median = np.median(lengths)
    std = np.std(lengths)
    if std > median:
        std = median
    else:
        pass
    for item in lengths:
        if abs(item - median) > (2. * std):
            lengths.remove(item)
        else:
            pass
    return int(np.mean(lengths)), int(np.std(lengths)), int(np.mean(rd))


def calc_cov(bam_name, start, stop):
    """
    calculates average coverage
    """
    bam = pysam.AlignmentFile(bam_name)
    # get chromosome names
    nms = []
    for i in bam.header['SQ']:
        if 'scaffold' in i['SN']:
            pass
        else:
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
    ustream = 0
    dstream = 0
    if chrom not in chrom_sizes.keys():
        raise Exception('Chromosome names do not match TE annotation')
    else:
        pass

    if (start - 2000) > 0:
        ustart = (start - 2000)
    else:
        ustart = 0

    if (stop + 2000) < chrom_sizes[chrom]:
        dstop = stop + 2000
    else:
        dstop = chrom_sizes[chrom]

    ul = start - ustart
    dl = dstop - stop
    l = stop - start

    for read in bam.pileup(chrom, start, stop):
        te += read.n
    for read in bam.pileup(chrom, ustart, start):
        ustream += read.n
    for read in bam.pileup(chrom, stop, dstop):
        dstream += read.n
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


def check_bam(bam, p, make_new_index=False):
    """
    Sort and index bam file
    returns dictionary of chromosome names and lengths
    """
    # check if sorted
    test_head = pysam.AlignmentFile(bam, 'rb')
    chrom_sizes = {}
    p = str(p)
    for i in test_head.header['SQ']:
        chrom_sizes[i['SN']] = int(i['LN'])
    try:
        test_head.header['HD']['SO']
    except KeyError:
        print '  sorting bam file'
        pysam.sort('-@', p, bam, 'sorted.temp')
        os.remove(bam)
        os.rename('sorted.temp.bam', bam)
    else:
        if test_head.header['HD']['SO'] == 'coordinate':
            pass
        else:
            print '  sorting bam file'
            pysam.sort('-@', p, bam, 'sorted.temp')
            os.remove(bam)
            os.rename('sorted.temp.bam', bam)
    test_head.close()
    # check if indexed
    if '{}.bai'.format(bam) in os.listdir('.') and make_new_index is False:
        pass
    else:
        print '  indexing bam file'
        pysam.index(bam)
    return chrom_sizes


def check_multi_te_deletion(coords, te_file):
    """
    If region spanned by split or disc reads is much larger than te,
    check if there are multiple TEs in the region that could all be deleted
    """
    coords = [str(x) for x in coords]
    tes = te_file.intersect(pybedtools.BedTool(" ".join(coords), from_string=True), nonamecheck=True)
    te_start_stop = []
    for i in tes:
        start = int(i[1])
        stop = int(i[2])
        te_start_stop.append([start, stop])
    merged = merge_intervals(te_start_stop)
    length = 0
    for i in merged:
        length += i[1]-i[0]
    return length


def merge_intervals(coords):
    """[[start,stop]...[start,stop]]"""
    copy = sorted(list(coords))
    for x in range(len(copy)-1):
        if len(copy) <= x+1:  # list shortened in recursion
            break
        if _overlap(copy[x][0], copy[x][1], copy[x+1][0], copy[x+1][1], 0) is True:
            start = min(copy[x][0], copy[x+1][0])
            stop = max(copy[x][1], copy[x+1][1])
            copy[x] = [start, stop]  # update copy
            del copy[x+1]  #remove merged
            copy = merge_intervals(copy)
        else:
            pass
    return copy


def determine_overlaps(coords, te_file, te_length, overlap, gapsize, target_te_overlap, target_gap_span, del_counts):
    """
    check if read spans enough of TE to call a deletion
    check that enough of TE spans gapsize
    """
    te_overlap = overlap / te_length
    read_overlap = te_length / gapsize
    if te_overlap < target_te_overlap:  # not enough TE covered to call deletion
        return False
    if read_overlap > target_gap_span:
        return True
    elif del_counts[",".join(map(str, coords))] > 1:
        multi_len = check_multi_te_deletion(coords, te_file)
        if (multi_len / gapsize) > target_gap_span:
            return True
        else:
            return False
    else:
        return False


def annotate_deletions(inp, acc, num_reads, bam, mn, p, te_file, del_counts):
    """
    Calls deletions where the gap between paired reads is at
    least 80 percent the length of the TE
    and there are either:
       num_reads split/disc read spanning the TE
    or
       coverage at TE is 1/10 the coverage in surrounding area and there are num_reads/2 reads
    """
    x = 0
    tes = {}
    written_tes = []
    chrom_sizes = check_bam(bam, p)
    allreads = pysam.AlignmentFile(bam, 'rb')

    with open(inp, 'r') as infile, open('deletions_{}.bed'.format(acc), 'w+') as outfile, open('deletion_reads_{}.txt'.format(acc), 'w+') as deletions_reads:
        for line in infile:
            line = line.rsplit()
            coords = [line[0], int(line[1]), int(line[2])]  # chr, start, stop
            te = [line[5], line[6], line[7], line[8], line[9]]  # chr, start, stop, strand, name
            name = te[4]
            length = int(te[2]) - int(te[1])
            overlap = int(line[12])
            gapsize = coords[2] - coords[1]
            read_type = line[4]
            read_name = line[3]
            if (gapsize <= 0) or (name in written_tes) or ((length-mn) > gapsize):
                pass
            else:
                if determine_overlaps(coords, te_file, length, overlap, gapsize, 0.8, 0.8, del_counts) is True:
                    if name not in tes.keys():
                        cov = get_coverages(coords[0], coords[1], coords[2], allreads, chrom_sizes)
                        tes[name] = [cov, 0, 0, [read_name]]  # coverage, split, disc, read_name (list)
                    else:
                        pass
                    if read_type == 'split':
                        tes[name][1] += 1
                        tes[name][3].append(read_name)
                    elif read_type == 'disc':
                        tes[name][2] += 1
                        tes[name][3].append(read_name)
                    else:
                        raise Exception('Incorrect read type information')
                    split = tes[name][1]
                    disc = tes[name][2]
                    total_reads = split + disc
                    if (tes[name][0] <= 0.1 and total_reads >= num_reads/2) or (total_reads >= num_reads):
                        data = (str(a) for a in te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=str(x)))
                        deletions_reads.write(">" + str(x) + "\t" + ",".join(tes[name][3]) + "\n")
                        x += 1
                        written_tes.append(name)
                    else:
                        pass
                else:
                    pass
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


def reorder_intersections(feature, read_count):
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
    disc_reads = int(feature[8])
    split_reads = int(feature[-2])
    total_reads = disc_reads + split_reads
    if total_reads >= int(read_count):
        if feature[3] == feature[13] and _overlap(int(feature[4]), int(feature[5]), int(feature[14]), int(feature[15]), 10) is True:
            feature = [chrom, start, stop, techrom, testart, testop, ','.join(reads), ','.join(names)]
            return feature
        else:
            pass
    else:
        pass


def check_name_sorted(bam, p):
    """
    Sort bam file by name if position sorted
    """
    test_head = pysam.AlignmentFile(bam, 'rb')
    p = str(p)
    try:
        test_head.header['HD']['SO']
    except KeyError:
        print '  sorting bam file'
        pysam.sort('-@', p, '-n', bam, 'sorted.temp')
        os.remove(bam)
        os.rename('sorted.temp.bam', bam)
        pysam.index(bam)
    else:
        if test_head.header['HD']['SO'] == 'queryname':
            pass
        else:
            print '  sorting bam file'
            pysam.sort('-@', p, '-n', bam, 'sorted.temp')
            os.remove(bam)
            os.rename('sorted.temp.bam', bam)
            pysam.index(bam)
    test_head.close()


def discover(options):
    """
    Discover TE insertions and deletions using read mapping information and TE annotation
    input concordant reads bam file must be position sorted
    imput split reads bam file must be name sorted
    TE annotation can be gzipped
    """
    print "Processing "+options.name
    print 'Estimating mean insert size and coverage'
    mn, std, rd_len = calc_mean(options.conc, options.proc)
    cov = calc_cov(options.conc, 100000, 120000)
    if cov <= 10:
        print '  Warning: coverage may not be sufficiently high to reliably discover polymorphic TE insertions'
    else:
        pass
    max_dist = (4*std) + mn
    print '\tmean insert size = {ins} bp, standard deviation = {std} bp\n\tcoverage = {cov}x\n\tread length = {rd} bp'.format(
        ins=mn, std=std, cov=cov, rd=rd_len)
    with open("tepy_discover_log_{}.txt".format(options.name), 'w+') as logfile:
        logfile.write('''Sample {sample}\nStart time {time}\nUsing TE annotation at {path}\nmean insert size = {ins} bp, standard deviation = {std} bp\ncoverage = {cov}x\nread length = {rd} bp\n'''.format(
                sample=options.name, time=ctime(),
                path=options.te, ins=mn, std=std, cov=cov, rd=rd_len))

    mask_chroms = options.mask.split(',')

    if options.strict is True:
        deletion_reads = int(cov/5) if (int(cov/5) > 10) else 10
        insertion_reads_low = int(cov/5) if (int(cov/5) > 10) else 10
        insertion_reads_high = int(cov/2) if (int(cov/2) > 10) else 10
        quality_filter_ins = 10
        quality_filter_dels = 10
    else:
        deletion_reads = int(cov/5) if (int(cov/5) > 4) else 4
        insertion_reads_low = int(cov/10) if (int(cov/10) > 2) else 2
        insertion_reads_high = int(cov/5) if (int(cov/5) > 2) else 2
        quality_filter_ins = 5
        quality_filter_dels = 0

    print 'Processing split reads'
    check_name_sorted(options.split, options.proc)
    split_unfiltered = pybedtools.BedTool(options.split).bam_to_bed().saveas().filter(lambda x: x[0] not in mask_chroms).saveas()
    split_unfiltered.filter(lambda x: int(x[4]) >= quality_filter_ins).saveas('split_hq.temp')
    split_unfiltered.filter(lambda x: int(x[4]) >= quality_filter_dels).saveas('split.temp')
    convert_split_pairbed('split_hq.temp', 'split_hq_bedpe.temp')
    convert_split_pairbed('split.temp', 'split_bedpe.temp')
    split_bedpe = pybedtools.BedTool('split_hq_bedpe.temp').each(append_origin, word='split').saveas().sort()
    split_bedpe_dels = pybedtools.BedTool('split_bedpe.temp').each(append_origin, word='split').saveas().sort()
    split_ins = split_bedpe.filter(lambda x: (abs(int(x[1]) - int(x[4])) > 5000) or (x[0] != x[3])).saveas()

    if options.discordant is False:
        print 'Finding discordant reads'
        filter_discordant(options.conc, max_dist, 'disc_bam.temp')
        pysam.sort('-@', str(options.proc), '-n', 'disc_bam.temp', 'disc_sorted')
        disc = pybedtools.BedTool('disc_sorted.bam')\
        .bam_to_bed(bedpe=True, mate1=True)\
        .filter(lambda x: x[0] not in mask_chroms).saveas()\
        .each(append_origin, word='disc').saveas()
    else:
        check_bam(options.discordant, options.proc)
        disc = pybedtools.BedTool(options.discordant)\
        .bam_to_bed(bedpe=True, mate1=True)\
        .filter(lambda x: x[0] not in mask_chroms).saveas()\
        .each(append_origin, word='disc').saveas()
    disc_split_dels = split_bedpe_dels.cat(disc, postmerge=False).sort().saveas('disc_split_dels.temp')
    disc_split_ins = split_ins.cat(disc, postmerge=False).sort().saveas('disc_split_ins.temp')

    print 'Processing TE annotation'
    te = pybedtools.BedTool(options.te).sort()
    disc_split_ins.pair_to_bed(te, f=0.80).saveas('intersect_ins.temp')

    if options.insertions is True:  # finding insertions only, so skip deletions
        pass
    else:
        print 'Finding deletions'
        create_deletion_coords(disc_split_dels, 'del_coords.temp')
        dels = pybedtools.BedTool('del_coords.temp').sort()
        dels.intersect(te, wo=True, sorted=True, nonamecheck=True).sort().saveas('deletions.temp')
        # need to intersect with TEs and record number of intersections for each coord to speed up processing later
        dels.intersect(te, c=True, sorted=True, nonamecheck=True).saveas('del_counts.temp')
        del_counts = {}
        with open('del_counts.temp', 'r') as infile:
            for line in infile:
                line = line.rsplit()
                crds = ",".join(line[:3])
                del_counts[crds] = int(line[-1])
        annotate_deletions('deletions.temp', options.name, deletion_reads, options.conc, mn, str(options.proc), te, del_counts)

    if options.deletions is True:  # finding deletions only, so skip insertions
        pass
    else:
        print 'Finding insertions'
        reorder('intersect_ins.temp', 'reorder_split.temp', 'forward_disc.temp', 'reverse_disc.temp')

        def merge_bed(infile, outfile):
            pybedtools.BedTool(infile).sort().merge(c='4,5,6,9,10,11,12,13,14',
                                                  o='collapse,collapse,collapse,distinct,collapse,count,collapse,collapse,collapse')\
            .saveas(outfile)

        file_pairs = [['reorder_split.temp','split_merged.temp'],
                      ['forward_disc.temp', 'forward_merged.temp'],
                      ['reverse_disc.temp', 'reverse_merged.temp']]

        for x in xrange(3):
            merge_bed(file_pairs[x][0], file_pairs[x][1])

        info = [['split_merged.temp', 'split_processed.temp', 'split'],
                ['forward_merged.temp', 'forward_processed.temp', 'disc_forward'],
                ['reverse_merged.temp', 'reverse_processed.temp', 'disc_reverse']]

        for x in xrange(3):
            process_merged(info[x][0], info[x][1], info[x][2])

        pybedtools.BedTool('forward_processed.temp').cat('reverse_processed.temp', postmerge=True,
            c='4,5,6,7,8,9,10',
            o='collapse,collapse,collapse,distinct,distinct,sum,distinct',
            d='200').sort().saveas('condensed_disc.temp')

        process_merged_disc('condensed_disc.temp', 'processed_disc.temp', 2, (mn+std), rd_len)
        pybedtools.BedTool('split_processed.temp').filter(lambda x: insertion_reads_high <= int(x[8])).saveas().each(lambda x: x[:-2]).moveto('high.temp')
        disc_split = pybedtools.BedTool('split_processed.temp')\
        .filter(lambda x: insertion_reads_high > int(x[8]))\
        .saveas().sort()\
        .intersect('processed_disc.temp', wo=True, nonamecheck=True)\
        .each(reorder_intersections, read_count=insertion_reads_low).saveas().sort().saveas()
        if len(disc_split) > 0:
            disc_split.cat('high.temp', postmerge=False).saveas().sort().saveas().moveto('insertions.temp')
            nm = 'insertions.temp'
        else:
            nm = 'high.temp'
        separate_reads(nm, 'insertions_{}.bed'.format(options.name), 'insertion_reads_{}.txt'.format(options.name))
    with open("tepy_discover_log_{}.txt".format(options.name), 'a') as logfile:
        logfile.write("tepy-discover finished normally at {}".format(ctime()))

    if options.keep is False:
        temp = glob('./*.temp')
        for i in temp:
            os.remove(i)
        os.remove('disc_sorted.bam')
