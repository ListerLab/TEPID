# usage:
# $ python annotate_ins.py a <accession_name> f <feature>
# where feature is gene or TE
# Finds insertion sites where there are reads at each end, refines insertion coordinates based on position of all reads

from sys import argv


def checkArgs(arg1, arg2):
    """
    arg1 is short arg, eg h
    arg2 is long arg, eg host
    """
    args = argv[1:]
    if arg1 in args:
        index = args.index(arg1)+1
        variable = args[index]
        return variable
    elif arg2 in args:
        index = args.index(arg2)+1
        variable = args[index]
        return variable
    else:
        variable = raw_input("\nEnter {arg2}: ".format(arg2=arg2))
        return variable


def annotate(collapse_file, insertion_file, id_file, accession_name):  # problem here
    """
    Find insertion coordinates and TE orientation. Adds unique ID: <accession_name>_<number>
    """
    with open(collapse_file, 'r') as infile, open(insertion_file, 'w+') as outfile, open(id_file, 'w+') as unique_id_file:
        outfile.write('ins_chr\tins_start\tins_stop\tins_strand\tAGI\tID\tref_chr\tref_start\tref_stop\tref_strand\n')
        i, lines = get_len(infile)
        ident = 0
        for x in range(i):
            line = lines[x]
            line = line.rsplit()
            chrom = line[0]
            start = line[1]
            stop = line[2]
            strand = line[3]
            te_name = line[4]
            orientation = line[5]
            mate = line[12]
            mate = mate.split(',')
            te_reads = line[11]
            te_reads = te_reads.split(',')
            reference = [line[7], line[8], line[9], line[10]]  # reference chrom, start, stop, strand
            pair = find_next(lines, i, x, int(chrom), strand, int(start), int(stop), te_name)
            if pair is False:
                pass  # no reads at opposite end, do not include in annotation
            else:
                pair_start = pair[0]
                pair_mates = pair[1]
                next_read_names = pair[2]
                mate = pair_mates + mate
                te_reads = next_read_names + te_reads
                outfile.write('{chr}\t{start}\t{stop}\t{orient}\t{name}\t{id}\t{ref}\n'.format(chr=chrom,
                                                                                               start=stop,
                                                                                               stop=pair_start,
                                                                                               orient=orientation,
                                                                                               name=te_name,
                                                                                               id=accession_name + '_' + str(ident),
                                                                                               ref='\t'.join(reference)))
                unique_id_file.write('>{id},{te},{reads},{mates}\n'.format(id=accession_name + '_' + str(ident),
                                                                           te=te_name,
                                                                           reads='|'.join(te_reads),
                                                                           mates='|'.join(mate)))  # Can add consensus sequence later
                ident += 1


# As file is processed top to bottom, sorted by coords, + will come up first. This will avoid identifying each insertion twice (once for each end)
def find_next(lines, i, x, chrom, strand, start, stop, te_name):
    """
    Find next read linked to same TE. Looks in 50 bp window.
    """
    while True:
        line = lines[x+1]
        line = line.rsplit()
        next_chrom = int(line[0])
        next_start = int(line[1])
        next_stop = int(line[2])
        next_strand = line[3]
        next_te_name = line[4]
        next_mate = line[11]
        next_mate = next_mate.split(',')
        next_te_reads = line[12]
        next_te_reads = next_te_reads.split(',')
        if strand != next_strand and te_name == next_te_name and chrom == next_chrom and stop <= next_start:
            return next_start, next_te_reads, next_mate
        elif stop + 50 < next_start:
            return False
        else:
            x += 1
            if x >= i:
                return False

accession = checkArgs('a', 'accession')
f = checkArgs('f', 'feature')

annotate('merged_{feat}_{acc}.bed'.format(acc=accession, feat=f),
         'insertions_{feat}_{acc}.bed'.format(acc=accession, feat=f),
         'id_{feat}_{acc}.fa'.format(acc=accession, feat=f),
         accession)
