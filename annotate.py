from subprocess import call


"""
1. Collapse overlapping coordinates mapping to the same TE
    - keep names of reads making up collapsed coordinates
2. define insertion point coordinates
    - Reads map to same te, different ends
    - reads are on different strands
    - reads are close together
3. Define orientation (strand)
    - which read maps to which end of TE
4. Write new file with:
    - insertion coordinates
    - orientation of TE
    - TE name
    - unique ID:
        - references read names making up annotation
        - can be used to compare between accessions
"""


def overlap(start1, stop1, start2, stop2):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    for y in xrange(start2, stop2):
        if start1 <= y <= stop1:
            return True
        else:
            pass


def get_len(infile):
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    return i, lines


def reorder(insert_file, reordered_file):  # should combine this with merge
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
            remaining = [field[6], field[13], field[14], field[15]]  # read name, TE strand, TE name, TE family
            te_coords = {'chrom': field[10], 'start': field[11], 'stop': field[12], 'strand': field[13], 'name': field[14]}
            if overlap(int(read1['start']), int(read1['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                te_read = read1
                dna_read = read2
            elif overlap(int(read2['start']), int(read2['stop']), int(te_coords['start']), int(te_coords['stop'])) is True:
                te_read = read2
                dna_read = read1
            else:
                raise Exception('check coords')
            merged = merge(x, lines, dna_read, remaining, outfile)  # need to sort first


def merge(x, lines, start, stop, te_name, strand, outfile):
    """
    merge coordinates that overlap where TE name and strand is same
    takes position in file (x), file (lines), start and stop coordinates, strand and te name, outfile to write merged lines.
    return list of read names, merged coordinates, and new position in file to start from
    """
    with open(sorted_file, 'r') as infile:
        mark = False
        i, lines = get_len(infile)
        x = 0
        while True:
            line = lines[x]
            field = line.rsplit()
            start = field[1]
            end = field[2]
            chrom = field[0]
            strand = field[3]
            te_name = field[10]
            nextline = lines[x+1]
            nextfield = nextline.rsplit()
            next_dna_start = nextfield[1]
            next_dna_end = nextfield[2]
            next_dna_chrom = nextfield[0]
            next_dna_strand = nextfield[3]
            next_te_start = nextfield[5]
            next_te_end = nextfield[6]
            next_te_chrom = nextfield[4]
            next_te_strand = nextfield[7]
            next_te_name = nextfield[10]
            if strand == next_dna_strand and te_name == next_te_name and chrom = next_dna_chrom:
                if overlap(next_dna_start, next_dna_end, start, stop) is True:
                    merge coordinates
                    x += 1
                else:
                    # write merged data to file, then start again from mark
                    # outfile.write('{chr1}\t{start1}\t{stop1}\t{strand1}\t{chr2}\t{start2}\t{stop2}\t{strand2}\t{remain}\n'.format(chr1=dna_read['chrom'],
                    #                                                                                                               start1=dna_read['start'],
                    #                                                                                                               stop1=dna_read['stop'],
                    #                                                                                                               strand1=dna_read['strand'],
                    #                                                                                                               chr2=te_read['chrom'],
                    #                                                                                                               start2=te_read['start'],
                    #                                                                                                               stop2=te_read['stop'],
                    #                                                                                                               strand2=te_read['strand'],
                    #                                                                                                               remain='\t'.join(remaining)))
                    marked_line = lines[mark+1]
                    marked_field = line.rsplit()
                    chrom = marked_field[0]
                    strand = marked_field[3]
                    te_name = marked_field[10]
                    x = mark
                    mark = False
            elif mark is False:
                mark = x
                x += 1
                if x >= i:
                    break
            else:
                x += 1
                if x >= i:
                    break


def annotate(collapse_dict):
    """
    Find insertion coordinates and TE orientation.
    """
    with open(collapse_file, 'r') as infile:
        lines = []
        for i, l in enumerate(infile):
            lines.append(l)
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            dna_read = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
            te_read = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
            te_coords = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[], 'name': line[]}
            pair = findNext(lines, x, dna_read['chrom'], dna_read['start'], dna_read['stop'], te_read['name'])
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                pair_dna = pair[0]
                pair_te = pair[1]
                insertion = {'chrom': pair_dna['chrom'], 'start': dna_read['stop'], 'stop': pair_dna['start'], 'strand': te_read['strand']}
                outfile.write('{chr}\t{start}\t{stop}\t{strand}\t{name}\t{id}\n'.format(chr=insertion['chrom'],
                                                                                        start=insertion['start'],
                                                                                        stop=insertion['stop'],
                                                                                        strand=insertion['strand'],
                                                                                        name=te_coords['name'],
                                                                                        id=))  # need to generate unique id
                else:
                    pass


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
def findNext(lines, x, chrom, strand, start, stop, te):
    """
    Find next read linked to same TE.
    """
    while True:
        line = lines[x+1]
        dna_read = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
        te_read = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
        te_name = line[]
        if strand == '+' and stop < dna_read['start'] + 200:
            return False
        elif strand == '-' and start > dna_read['stop'] - 200:
            return False
        elif strand != dna_read['strand'] and te == te_read['name']:
            return dna_read, te_read
        else:
            x += 1


reorder('{b}.bed'.format(b=accession_name, 'ordered_{b}.bed'.format(b=accession_name)))
call('sort -n1,1 -nk2,2 ordered_{b}.bed > sorted_{b}.bed'.format(b=accession_name), shell=True)
annotate('sorted_{b}.bed'.format(b=accession_name))  # insert filename
