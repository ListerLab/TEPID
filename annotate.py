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


def collapse(insert_file):
    with open(insert_file, 'r') as infile, open('insertions.tsv', 'w+') as outfile:
        lines = []
        for i, l in enumerate(infile):
            lines.append(l)
        for x in range(i):
            line = lines[x]
            field = line.rsplit()
            chrom1 = 
            start1 = 
            stop1 = 
            strand1 = 
            chrom2 = 
            start2 = 
            stop2 = 
            strand2 = 
            te = 
            te_chrom = 
            te_start = 
            te_stop = 
            te_strand = 


def annotate(collapse_file):
    with open(collapse_file, 'r') as infile:
        lines = []
        for i, l in enumerate(infile):
            lines.append(l)
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            read1 = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
            read1 = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
            te_coords = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[], 'name': line[]}
            # Find out which read is in TE
            if te_coords['start'] < read1['start'] > te_coords['stop'] and read1['chrom'] == te_coords['chrom']:
                te_read = read1
                dna_read = read2
            elif te_coords['start'] < read2['start'] > te_coords['stop'] and read2['chrom'] == te_coords['chrom']:
                te_read = read2
                dna_read = read1
            else:
                raise Exception('check coords')
            # Define TE orientation
            if dna_read['strand'] == te_read['strand']:
                orientation = te_read['strand']
            else:
                if te_read['strand'] == '-':
                    orientation = '+'
                else:
                    orientation = '-'
            # Find other end of TE insertion point
            pair = findNext(lines, x, dna_read['chrom'], dna_read['start'], dna_read['stop'], te_read['name'])
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                pair_dna = pair[0]
                pair_te = pair[1]
                if dna_read['strand'] == '+':  # dna_read is first
                    insertion = {'chrom': dna_read['chrom'], 'start': dna_read['stop'], 'stop': pair_dna['start'], 'strand': orientation}
            outfile.write('{chr}\t{start}\t{stop}\t{strand}\t{name}\t{id}\n'.format(chr=insertion['chrom'],
                                                                                    start=insertion['start'],
                                                                                    stop=insertion['stop'],
                                                                                    strand=insertion['strand'],
                                                                                    name=te_coords['name'],
                                                                                    id=))


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
# might cause a problem when dna read is in second position
# Could reformat file so dna read is always in first position and sort by chrom, start
def findNext(lines, x, chrom, strand, start, stop, te):
    while True:
        line = lines[x+1]
        read1 = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
        read1 = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[]}
        te_coords = {'chrom': line[], 'start': line[], 'stop': line[], 'strand': line[], 'name': line[]}
        if strand == '+' and stop < dna_read['start'] + 200:
            return False
        elif strand == '-' and start > dna_read['stop'] - 200:
            return False
        elif strand is not dna_read['strand'] and te == te_read['name']:
            return dna_read, te_read
        else:
            x += 1


@collapse
annotate(bedfile)
