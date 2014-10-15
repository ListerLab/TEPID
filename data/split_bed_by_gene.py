
def main(bedfile, splitcol, out_prefix):
    bedfh = open(bedfile)
    outfiles = {}
    for line in bedfh:
        # turn line into a list of strings
        line = line.strip().split('\t')
        # grab out the gene
        gene = line[splitcol]
        try:
            # try grabbing this gene's file from the dict
            ofh = outfiles[gene]
        except KeyError:
            # open file and store in dict if we don't
            ofh = open(out_prefix + gene, "w")
            outfiles[gene] = ofh
        # write the line, inluding newline we stripped, to the gene's file
        ofh.write('\t'.join(line) + '\n')
    # close the files
    for fh in outfiles.values():
        fh.close()
    bedfh.close()

if __name__ == "__main__":
    import sys
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
