# usage:
# python inverse_del.py <names> <merged_deletions> <outfile> <reference_sample>
# requires list of all sample names


def filter_del(inf, master, outf, ref):
    """
    Take bedfile containing all TE deletions and create
    polymorphic TE file with coordinates of TE and
    list of accessions that contain the TE
    """
    with open(inf, 'r') as infile, open(outf, 'w+') as outfile:
        for line in infile:
            line = line.rsplit()
            accessions = line[6]
            accessions = accessions.split(',')
            coords = line[:5]
            temp = [ref]
            for item in master:
                if item not in accessions:
                    temp.append(item)
                else:
                    pass
            info = '\t'.join(coords) + '\t' + ','.join(temp) + '\n'
            outfile.write(info)


def create_names_list(inf):
    names = []
    with open(inf, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            names.append(line[0])
    return names

if __name__ == "__main__":
    import sys
    master = create_names_list(sys.argv[1])
    filter_del(sys.argv[2], master, sys.argv[3], sys.argv[4])
