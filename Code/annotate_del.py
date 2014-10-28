from __future__ import division
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


def annotate_deletions(accession):
    """
    Calls deletions where the gap between paired reads is at
    least 75 percent the length of the TE
    """
    with open("{acc}_deletions_temp.bed".format(acc=accession), 'r') as infile, open("deletions_{acc}.bed".format(acc=accession), 'w+') as outfile:
        x = 0
        for line in infile:
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
                    ident = 'del_{acc}_{x}'.format(acc=accession, x=x)
                    data = map(str, te)
                    outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                    x += 1
                    name = te[4]
                else:
                    if name != te[4]:
                        ident = 'del_{acc}_{x}'.format(acc=accession, x=x)
                        data = map(str, te)
                        outfile.write('{te}\t{id}\n'.format(te='\t'.join(data), id=ident))
                        x += 1
                        name = te[4]
                    else:
                        name = te[4]

a = checkArgs('a', 'accession')

annotate_deletions(a)
