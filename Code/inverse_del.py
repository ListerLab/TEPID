
def filter_del(inf, master):
    """
    Take bedfile containing all TE deletions and create
    polymorphic TE file with coordinates of TE and
    list of accessions that contain the TE
    """
    with open(inf, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            accessions = line[6]
            accessions = accessions.split(',')
            coords = line[:5]
            ins = []
            for item in accessions:
                if item not in master:
                    ins.append(item)
                    ins.append('Col-0')  # will be in all because reference
                else:
                    pass
            info = '\t'.join(coords) + '\t' + ','.join(ins) + '\n'
            print info

