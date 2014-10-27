# usage:
# python separate_read_info.py <acc>
# where <acc> is accession name


def splitfile(acc):  # this needs work
    with open('insertions_TE_{a}_split_reads.bed'.format(a=acc), 'r') as infile, open('single_break_temp.bed', 'w+') as single, open('double_break_temp.bed', 'w+') as double:
        for line in infile:
            field = line.rsplit()
            count = int(field[11])
            data = field[:11]
            if count == 1:
                single.write("{data}\n".format(data='\t'.join(data)))
            elif count == 2:
                double.write("{data}\n".format(data='\t'.join(data)))
            else:
                pass

if __name__ == "__main__":
    import sys
    splitfile(sys.argv[1])
