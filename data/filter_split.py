
def filter_lines(acc):
    with open('{a}_merged.split.bed'.format(a=acc), 'r') as infile, open('{a}_filtered.split.bed'.format(a=acc), 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            chrom = int(line[0])
            start = int(line[1])
            stop = int(line[2])
            count_start = int(line[3])
            count_stop = int(line[4])
            if count_stop == 1 and count_start == 1:
                pass
            elif count_start == 1:
                outfile.write('{chr}\t{brk}\t{brk}\n'.format(chr=chrom, brk=start))
            elif count_stop == 1:
                outfile.write('{chr}\t{brk}\t{brk}\n'.format(chr=chrom, brk=stop))
            else:
                pass

if __name__ == '__main__':
    import sys
    filter_lines(sys.argv[1])
