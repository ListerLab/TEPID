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
            remaining = [field[6], field[13], field[14], field[15]]  # read name, TE strand, TE name, TE family
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
            outfile.write('{chr1}\t{start1}\t{stop1}\t{strand1}\t{chr2}\t{start2}\t{stop2}\t{strand2}\t{remain}\t{mate}\n'.format(chr1=dna_read['chrom'],
                                                                                                                                  start1=dna_read['start'],
                                                                                                                                  stop1=dna_read['stop'],
                                                                                                                                  strand1=dna_read['strand'],
                                                                                                                                  chr2=te_read['chrom'],
                                                                                                                                  start2=te_read['start'],
                                                                                                                                  stop2=te_read['stop'],
                                                                                                                                  strand2=te_read['strand'],
                                                                                                                                  remain='\t'.join(remaining),
                                                                                                                                  mate=mate))


def merge(sorted_file, output_file):
    """
    merge coordinates that overlap where TE name and strand is same
    takes sorted output from reorder function and name of output file.
    return file with of read names, merged coordinates, and new position in file to start from
    """
    with open(sorted_file, 'r') as infile, open(output_file, 'w+') as outfile:
        i, lines = get_len(infile)
        x = 1
        readlist = []
        done = True
        while True:
            if done is True:
                line = lines[x]
                field = line.rsplit()
                start = field[1]
                end = field[2]
                chrom = field[0]
                strand = field[3]
                te_name = field[10]
                readname = field[8]
                te_strand = field[7]
                mate = field[12]  # needs to be a list matched with read names
            else:
                pass
            nextline = lines[x+1]
            nextfield = nextline.rsplit()
            next_dna_start = nextfield[1]
            next_dna_end = nextfield[2]
            next_dna_chrom = nextfield[0]
            next_dna_strand = nextfield[3]
            next_read_name = nextfield[8]
            next_te_name = nextfield[10]
            if strand == next_dna_strand and te_name == next_te_name and chrom == next_dna_chrom:
                if overlap(int(next_dna_start), int(next_dna_end), int(start), int(end)) is True:
                    readlist.append(next_read_name)
                    if start > next_dna_start:
                        start = next_dna_strand
                    elif end < next_dna_end:
                        end = next_dna_end
                    else:
                        pass
                    x += 1
                    done = False
                else:
                    if readlist:
                        reads = ','.join(readlist)
                    else:
                        reads = readname
                    outfile.write('{chr}\t{start}\t{end}\t{strand}\t{te}\t{orient}\t{reads}\t{mate}\n'.format(chr=chrom,
                                                                                                              start=start,
                                                                                                              end=end,
                                                                                                              strand=strand,
                                                                                                              te=te_name,
                                                                                                              orient=te_strand,
                                                                                                              reads=reads,
                                                                                                              mate=mate))
                    done = True
                    x += 1
                    readlist = []
            else:
                if readlist:
                    reads = ','.join(readlist)
                else:
                    reads = readname
                outfile.write('{chr}\t{start}\t{end}\t{strand}\t{te}\t{orient}\t{reads}\t{mate}\n'.format(chr=chrom,
                                                                                                          start=start,
                                                                                                          end=end,
                                                                                                          strand=strand,
                                                                                                          te=te_name,
                                                                                                          orient=te_strand,
                                                                                                          reads=reads,
                                                                                                          mate=mate))
                done = True
                x += 1
                readlist = []
                if x >= i:
                    break


def annotate(collapse_dict, insertion_file, id_file):
    """
    Find insertion coordinates and TE orientation.
    """
    with open(collapse_file, 'r') as infile, open(insertion_file, 'w+') as outfile, open(id_file, 'w+') as unique_id_file:
        lines = []
        for i, l in enumerate(infile):
            lines.append(l)
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            chrom = line[0]
            start = line[1]
            stop = line[2]
            strand = line[3]
            te_name = line[4]
            reads = line[6]
            te_reads = line[7]
            pair = findNext(lines, x, chrom, start, stop, te_name)
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                outfile.write('{chr}\t{start}\t{stop}\t{strand}\t{name}\t{id}\n'.format(chr=chrom,
                                                                                        start=stop,
                                                                                        stop=pair,
                                                                                        strand=strand,
                                                                                        name=te_name,
                                                                                        id=))  # need to generate unique id
                unique_id_file.write()  # write id, list of read names and which is the TE read, separated in fasta style. Can later add consensus sequence
                else:
                    pass


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
def findNext(lines, x, chrom, strand, start, stop, te_name):
    """
    Find next read linked to same TE. Looks in 200 bp window.
    """
    while True:
        line = lines[x+1]
        line = line.rsplit()
        next_chrom = line[0]
        next_start = line[1]
        next_stop = line[2]
        next_strand = line[3]
        next_te_name = line[4]
        next_reads = line[6]
        next_te_reads = line[7]
        if strand == '+' and stop < next_start + 200:
            return False
        elif strand == '-' and start > next_stop - 200:
            return False
        elif strand != next_strand and te_name == next_te_read:
            return next_start
        else:
            x += 1


reorder('{b}.bed'.format(b=accession_name, 'ordered_{b}.bed'.format(b=accession_name)))
call('sort -n1,1 -nk2,2 ordered_{b}.bed > sorted_{b}.bed'.format(b=accession_name), shell=True)
merge('sorted_{b}.bed'.format(b=accession_name), 'merged_{b}.bed'.format(b=accession_name))
annotate('merged_{b}.bed'.format(b=accession_name), 'insertions_{b}.bed'.format(b=accession_name))
