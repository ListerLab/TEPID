
def filter_lines(inp):
    """
    Takes human transposon annotation file from
    UCSC repeatmasker track and reformats as bedfile
    """
    for line in inp:
        line = line.rsplit()
        chrom = line[5]
        start = line[6]
        stop = line[7]
        strand = line[9]
        name = line[10]
        cls = line[11]
        fam = line[12]
        if 'Simple_repeat' in fam or 'Low_complexity' in fam or 'tRNA' in fam or 'rRNA' in fam:
            pass
        elif len(chrom) < 6:
            print ('{ch}\t{sta}\t{sto}\t{stra}\t{nm}\t{fam}\t{cls}'.format(ch=chrom,
                                                                           sta=start,
                                                                           sto=stop,
                                                                           stra=strand,
                                                                           nm=name,
                                                                           fam=fam,
                                                                           cls=cls))
        else:
            pass

if __name__ == "__main__":
    import sys
    filter_lines(sys.stdin)
