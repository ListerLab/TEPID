#! /usr/bin/env python

from argparse import ArgumentParser


parser = ArgumentParser(description='Genotype TE insertions')
group = parser.add_mutually_exclusive_group()
group.add_argument('-d', '--deletions', help='run on deletions', action='store_true', required=False, default=False)
group.add_argument('-i', '--insertions', help='run on insertions', action='store_true', required=False, default=False)
parser.add_argument('-a', '--ambiguous', help='ambiguous TE variants filename', required=True)
parser.add_argument('-m', '--merged', help='merged TE variants filename', required=True)
parser.add_argument('-s', '--samples', help='all sample names', required=True)
parser.add_argument('-r', '--reference', help='reference sample name', required=True)
options = parser.parse_args()

if options.deletions is True:
    val = 5  # number of columns
elif options.insertions is True:
    val = 7
else:
    raise Exception("Incorrect arguments")


def create_names_list(inf):
    names = []
    with open(inf, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            names.append(line[0])
    return names


def read_files_to_dict(f, val):
    d = {}
    with open(f, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            coords = '|'.join(line[:val])  # key common to ambiguous and merged insertions files
            d[coords] = line[-1].split(',')  # accession names
    return d


def invert_samples(samples, all_accessions, reference):
    i = [reference]
    for accession in all_accessions:
        if accession not in samples:
            i.append(accession)
    return i


def genotype(merged, ambiguous, all_accessions, reference):
    for key, value in merged.items():
        opposite_accessions = invert_samples(value, all_accessions, reference)
        try:
            ambiguous[key]
        except KeyError:
            ambiguous_accessions = []
        else:
            ambiguous_accessions = ambiguous[key]
        for i in opposite_accessions:
            if i in ambiguous_accessions:
                opposite_accessions.remove(i)
        data = key.split('|')
        te = data[-1]
        coords = data[:-1]
        print("\t".join(coords)+"\t"+te+"\t"+",".join(value)+"\t"+",".join(opposite_accessions))


accession = create_names_list(options.samples)
merged = read_files_to_dict(options.merged, val)
ambiguous = read_files_to_dict(options.ambiguous, val)
genotype(merged, ambiguous, accession, options.reference)