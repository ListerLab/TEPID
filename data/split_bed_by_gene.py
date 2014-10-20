
def main(bedfile, splitcol, out_prefix):
    """
    Takes input bedfile sorted by TE name and splits into
    different file for each TE name
    """
    with open(bedfile, 'r') as bedin:
        for line in bedin:
            field = line.rsplit()
            te = field[splitcol]
            try:
                prev_te
            except NameError:
                outfile = open("{pref}_{te}". format(pref=out_prefix, te=te), 'w')
                outfile.write(line)
                prev_te = te
            else:
                if te == prev_te:
                    outfile.write(line)
                    prev_te = te
                else:
                    outfile.close()
                    outfile = open("{pref}_{te}". format(pref=out_prefix, te=te), 'w')
                    outfile.write(line)
                    prev_te = te
        outfile.close()

if __name__ == "__main__":
    import sys
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
