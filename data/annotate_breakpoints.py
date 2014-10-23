
def get_len(infile):
    """returns number of lines in file and all lines as part of list"""
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    return i, lines


def annotate_single(acc):
    """adds breakpoint coordinates to insertion"""
    with open('single_break.bed', 'r') as infile, open('insertions_{a}_temp.bed'.format(a=acc), 'a+') as outfile:
        for line in infile:
            line = line.rsplit()
            coords = line[11:14]
            data = line[3:11]
            outfile.write('{coords}\t{data}\n'.format(coords='\t'.join(coords), data='\t'.join(data)))


def annotate_double(acc):
    """adds breakpoint coordinates to insertions with two breakpoints"""
    with open('double_break.bed', 'r') as infile, open('insertions_{a}_temp.bed'.format(a=acc), 'a+') as outfile:
        i, lines = get_len(infile)
        x = 0
        while x < i:
            line = lines[x].rsplit()
            chrom = int(line[11])
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

if __name__ == '__main__':
    import sys
    annotate_single(sys.argv[1])
    annotate_double(sys.argv[1])
