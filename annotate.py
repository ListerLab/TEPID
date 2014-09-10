from subprocess import call
from sys import argv

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
        matelist = []
        done = True
        while True:
            if done is True:
                line = lines[x]
                field = line.rsplit()
                start = field[1]
                end = field[2]
                chrom = field[0]
                strand = field[3]
                te_name = field[9]
                readname = field[11]
                te_strand = field[4]
                reference = [field[10], field[5], field[6], field[7], field[8]]
                mate = field[12]
            else:
                pass
            nextline = lines[x+1]
            nextfield = nextline.rsplit()
            next_dna_start = nextfield[1]
            next_dna_end = nextfield[2]
            next_dna_chrom = nextfield[0]
            next_dna_strand = nextfield[3]
            next_read_name = nextfield[11]
            next_te_name = nextfield[9]
            next_mate = nextfield[12]
            if strand == next_dna_strand and te_name == next_te_name and chrom == next_dna_chrom:
                if overlap(int(next_dna_start), int(next_dna_end), int(start), int(end)) is True:
                    readlist.append(next_read_name)
                    matelist.append(next_mate)
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
                        matelist.append(mate)
                        readlist.append(readname)
                        reads = ','.join(readlist)
                        mates = ','.join(matelist)
                    else:
                        reads = readname
                        mates = mate
                    outfile.write('{chr}\t{start}\t{end}\t{strand}\t{te}\t{orient}\t{ref}\t{reads}\t{mate}\n'.format(chr=chrom,
                                                                                                                     start=start,
                                                                                                                     end=end,
                                                                                                                     strand=strand,
                                                                                                                     te=te_name,
                                                                                                                     orient=te_strand,
                                                                                                                     ref='\t'.join(reference),
                                                                                                                     reads=reads,
                                                                                                                     mate=mates))
                    done = True
                    x += 1
                    readlist = []
                    matelist = []
            else:
                if readlist:
                    matelist.append(mate)
                    readlist.append(readname)
                    reads = ','.join(readlist)
                    mates = ','.join(matelist)
                else:
                    reads = readname
                    mates = mate
                outfile.write('{chr}\t{start}\t{end}\t{strand}\t{te}\t{orient}\t{ref}\t{reads}\t{mate}\n'.format(chr=chrom,
                                                                                                                 start=start,
                                                                                                                 end=end,
                                                                                                                 strand=strand,
                                                                                                                 te=te_name,
                                                                                                                 orient=te_strand,
                                                                                                                 ref='\t'.join(reference),
                                                                                                                 reads=reads,
                                                                                                                 mate=mates))
                done = True
                x += 1
                readlist = []
                matelist = []
                if x >= i:
                    break


def annotate(collapse_file, insertion_file, id_file):
    """
    Find insertion coordinates and TE orientation.
    """
    with open(collapse_file, 'r') as infile, open(insertion_file, 'w+') as outfile, open(id_file, 'w+') as unique_id_file:
        outfile.write('ins_chr\tins_start\tins_stop\tins_strand\tAGI\tID\tref_chr\tref_start\tref_stop\tref_strand\n')
        lines = []
        for i, l in enumerate(infile):
            lines.append(l)
            ident = 0
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            chrom = line[0]
            start = line[1]
            stop = line[2]
            strand = line[3]
            te_name = line[4]
            orientation = line[5]
            mate = line[12]
            mate = mate.split(',')
            te_reads = line[11]
            te_reads = te_reads.split(',')
            reference = [line[7], line[8], line[9], line[10]]  # reference chrom, start, stop, strand
            pair = find_next(lines, i, x, int(chrom), strand, int(start), int(stop), te_name)
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                pair_start = pair[0]
                pair_mates = pair[1]
                next_read_names = pair[2]
                mate = pair_mates + mate
                te_reads = next_read_names + te_reads
                outfile.write('{chr}\t{start}\t{stop}\t{orient}\t{name}\t{id}\t{ref}\n'.format(chr=chrom,
                                                                                               start=stop,
                                                                                               stop=pair_start,
                                                                                               orient=orientation,
                                                                                               name=te_name,
                                                                                               id=ident,
                                                                                               ref='\t'.join(reference)))
                unique_id_file.write('>{id},{te},{reads},{mates}\n'.format(id=ident,
                                                                           te=te_name,
                                                                           reads='|'.join(te_reads),
                                                                           mates='|'.join(mate)))  # Can add consensus sequence later
                ident += 1


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
def find_next(lines, i, x, chrom, strand, start, stop, te_name):
    """
    Find next read linked to same TE. Looks in 50 bp window.
    """
    while True:
        line = lines[x+1]
        line = line.rsplit()
        next_chrom = int(line[0])
        next_start = int(line[1])
        next_stop = int(line[2])
        next_strand = line[3]
        next_te_name = line[4]
        next_mate = line[11]
        next_mate = next_mate.split(',')
        next_te_reads = line[12]
        next_te_reads = next_te_reads.split(',')
        if strand != next_strand and te_name == next_te_name and chrom == next_chrom:
            return next_start, next_te_reads, next_mate
        elif stop + 50 < next_start:
            return False
        else:
            x += 1
            if x >= i:
                return False


accession_name = checkArgs('n', 'name')
reorder('{b}_TE_intersections.bed'.format(b=accession_name), 'intersections_ordered_{b}.bed'.format(b=accession_name))
call('sort -n1,1 -nk2,2 intersections_ordered_{b}.bed > intersections_sorted_{b}.bed'.format(b=accession_name), shell=True)
merge('intersections_sorted_{b}.bed'.format(b=accession_name), 'merged_{b}.bed'.format(b=accession_name))
annotate('merged_{b}.bed'.format(b=accession_name), 'insertions_{b}.bed'.format(b=accession_name), 'id_{b}.fa'.format(b=accession_name))
