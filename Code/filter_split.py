
def filter_lines(inp):
    for line in inp:
        line = line.rsplit()
        chrom = int(line[0])
        start = int(line[1])
        stop = int(line[2])
        count_start = int(line[3])
        count_stop = int(line[4])
        if count_stop == 1 and count_start == 1:
            pass
        elif count_start == 1:
            print ('{chr}\t{brk}\t{brk}\n'.format(chr=chrom, brk=start))
        elif count_stop == 1:
            print ('{chr}\t{brk}\t{brk}\n'.format(chr=chrom, brk=stop))
        else:
            pass

if __name__ == '__main__':
    import sys
    filter_lines(sys.stdin)
