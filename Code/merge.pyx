def merge_te_coords(infile, outfile, num_reads):
    """
    takes file containing reads coordinates
    that overlap annotated TEs and creates a 
    new file with merged coordinates
    only writes lines if there are at least
    num_split and num_disc split and disc reads
    """
    with open(infile, 'r') as inp:
        TE_dict = _create_te_dict(inp)
    with open(outfile, 'w+') as outf:
        for name in TE_dict.keys():
            _modify_coords(TE_dict[name])
            for key, value in TE_dict[name].items():
                chrom = value[0]
                start = value[1]
                stop = value[2]
                strand = value[5]
                te_coords = value[3]
                ref = '\t'.join(te_coords)
                reads = ','.join(value[4])
                mates = ','.join(value[6])
                read_count = len(value[7])
                breakpoint = value[8]
                if 'True' in breakpoint:
                    breakpoint = True
                else:
                    breakpoint = False
                if (read_count >= num_reads) or ((breakpoint == True )and (read_count > num_reads-10)):
                    outf.write('{ch}\t{sta}\t{sto}\t{name}\t{ref}\t{reads}\t{mates}\n'.format(ch=chrom,
                                                                                              sta=start,
                                                                                              sto=stop,
                                                                                              name=name,
                                                                                              ref=ref,
                                                                                              reads=reads,
                                                                                              mates=mates))
                else:
                    pass


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


def _modify_coords(inp):
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
                sd = [line[13]]
                breakpoints = [line[14]]
                _merge(chrom1, start1, stop1, inp, x, strands, reads, mates, ref_coords, skips, sd, breakpoints)
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
    sd = [line[13]]
    breakpoints = [line[14]]
    inp[x] = [chrom1, start1, stop1, ref_coords, reads, strands, mates, sd, breakpoints]


def _merge(chrom1, start1, stop1, d, x, strands, reads, mates, ref_coords, skips, sd, breakpoints):
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
                strands.append(line[3])
                mates.append(line[12])
                sd.append(line[13])
                breakpoints.append(line[14])
                del d[key]
                skips.append(key)
            else:
                start = start1
                stop = stop1
    try:
        start
    except NameError:  # all keys were skipped as they were merged with other reads already
        start = start1
        stop = stop1
        d[x] = [chrom1, start, stop, ref_coords, reads, strands, mates, sd, breakpoints]
        skips.append(x)
    else:
        strands = list(set(strands))
        d[x] = [chrom1, start, stop, ref_coords, reads, strands, mates, sd, breakpoints]
        skips.append(x)


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