from __future__ import division


def count_inserts(inf, outf, chrom):
    """
    Splits chromosome into 50 bins.
    Finds frequency of TE inserts in each bin.
    Returns bin start points and number of insertions in each bin as tsv file.
    """
    with open(inf, 'r') as infile, open(outf, 'w+') as outfile:
        chr_size = {1: 30427617, 2: 19698289, 3: 23459830, 4: 18585056, 5: 26975502}
        length = chr_size[chrom]
        bins = length / 50
        bins = int(bins)
        bins_dict = {}
        for x in range(50):
            bins_dict[x] = 0
        for line in infile:
            line = line.rsplit()
            ins_chr = line[0]
            ins_chr = int(ins_chr)
            ins_start = int(line[1])
            ins_end = line[2]
            if ins_chr == chrom:
                bin_no = which_bin(bins, 50, ins_start)
                bins_dict[bin_no] += 1
            else:
                pass
        for key, value in bins_dict.items():
            outfile.write('{v}\n'.format(v=value))
        # outfile.write('{keys}\n'.format(keys='\t'.join(map(str, bins_dict.keys()))))
        # outfile.write('{values}\n'.format(values='\t'.join(map(str, bins_dict.values()))))


def which_bin(size, number, inp):
    """
    Takes number of bins, bin size, and a number to evaluate.
    Finds which bin the input number goes into.
    Returns bin number.
    """
    if inp <= size:
        return 1
    else:
        for x in xrange(number):
            if (x*size) <= inp <= ((x+1)*size):
                return x
            else:
                pass

count_inserts('insertions_Aa-0.bed', 'chrom1.tsv', 1)
