#! /usr/bin/python

from argparse import ArgumentParser


parser = ArgumentParser(description='Invert the TE deletion calls to give a consistent data format between TE insertions and deletions')
parser.add_argument('-s', '--samples', help='list of all sample names', required=True)
parser.add_argument('-d', '--deletions', help='merged TEPID deletions', required=True)
parser.add_argument('-r', '--reference', help='reference sample name, eg Col-0', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)
options = parser.parse_args()


def filter_del(options):
    with open(options.deletions, 'r') as dels, open(options.output, 'w+') as outfile:
        master = [line.strip("\n") for line in open(options.samples, "r")]
        for line in dels:
            line = line.rsplit()
            accessions = line[5]
            accessions = accessions.split(',')
            coords = line[:4]
            temp = [options.reference]
            te = line[4]
            for item in master:
                if item not in accessions:
                    temp.append(item)
                else:
                    pass
            coords.pop(3)  # remove strand
            info = '\t'.join(coords) + '\t' + te + '\t' + ','.join(temp) + '\n'
            outfile.write(info)


if __name__ == "__main__":
    filter_del(options)
