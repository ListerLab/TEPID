from __future__ import division


def annotate_deletions(inp, acc):
    """
    Calls deletions where the gap between paired reads is at
    least 75 percent the length of the TE
    """
    x = 0
    for line in inp:
        line = line.rsplit()
        coords = [int(line[0]), int(line[1]), int(line[2])]  # chr, start, stop
        te = [line[4], line[5], line[6], line[7], line[8]]  # chr, start, stop, strand, name
        overlap = int(line[11])
        gapsize = coords[2] - coords[1]
        percentage = overlap / gapsize
        if percentage >= 0.75:
            try:
                name
            except NameError:
                ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                data = map(str, te)
                print ('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                x += 1
                name = te[4]
            else:
                if name != te[4]:
                    ident = 'del_{acc}_{x}'.format(acc=acc, x=x)
                    data = map(str, te)
                    print ('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                    x += 1
                    name = te[4]
                else:
                    name = te[4]
        else:
            pass


if __name__ == "__main__":
    import sys
    annotate_deletions(sys.stdin, sys.argv[1])
