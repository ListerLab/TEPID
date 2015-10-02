from __future__ import division
import os
import numpy as np
import pysam
import heapq


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


def calc_cov(bam_name, start, stop, p):
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


def check_bam(bam, p):
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
    if test_head.header['HD']['SO'] == 'coordinate':
        pass
    else:
        print '  sorting bam file'
        pysam.sort('-@', p, bam, 'sorted.temp')
        os.remove(bam)
        os.rename('sorted.temp.bam', bam)
    test_head.close()
    # check if indexed
    if '{}.bai'.format(bam) in os.listdir('.'):
        pass
    else:
        print '  indexing bam file'
        pysam.index(bam)
    return chrom_sizes


def annotate_deletions(inp, acc, num_reads, bam, mn, p):
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
                percentage = overlap / gapsize
                if percentage >= 0.8:
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
                    # 10% coverage and half the normally required reads
                    if (tes[name][0] <= 0.1 and total_reads >= num_reads/2) or (total_reads >= num_reads):
                        ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                        data = (str(x) for x in te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                        deletions_reads.write(">" + ident + "\t" + ",".join(tes[name][3]) + "\n")
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
    disc_reads = int(feature[8])
    split_reads = int(feature[-2])
    if disc_reads >= int(num_disc) or split_reads >= int(num_split):
        if feature[3] == feature[13] and _overlap(int(feature[4]), int(feature[5]), int(feature[14]), int(feature[15]), 10) is True:
            feature = [chrom, start, stop, techrom, testart, testop, ','.join(reads), ','.join(names)]
            return feature
        else:
            pass
    else:
        pass
