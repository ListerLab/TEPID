
def get_len(infile):
    """returns number of lines in file and all lines as part of list"""
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    return i, lines


def convert_split_pairbed(inp):
    i, lines = get_len(inp)
    x = 0
    while x < i:
        coords, read, strand = get_features(lines[x])
        x += 1
        next_coords, next_read, next_strand = get_features(lines[x])
        if next_read == read:
            mate = read.split('_')
            print ("{co}\t{nco}\t{read}\t{mt}\t{st1}\t{st2}".format(co='\t'.join(coords),
                                                                    nco='\t'.join(next_coords),
                                                                    read=mate[0],
                                                                    mt=mate[1],
                                                                    st1=strand,
                                                                    st2=next_strand))
            x += 1
        else:
            pass


def get_features(inp):
    line = inp.rsplit()
    coords = line[:3]
    read = line[3]
    strand = line[5]
    return coords, read, strand

if __name__ == "__main__":
    import sys
    convert_split_pairbed(sys.argv[1])
