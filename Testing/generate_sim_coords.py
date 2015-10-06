#! /usr/bin/env python

import random


class sim_data():

    def __init__(self):
        self.chromosomes = {}
        self.transposons = {}

    def read_genome_index(self, index):
        """
        reads and stores chromosome information from genome indexed by samtools faidx
        """
        with open(index, 'r') as infile:
            for line in infile:
                line = line.rsplit()
                if 'super' in line[0] or 'scaffold' in line[0]:
                    pass
                elif line[0].startswith('chr'):
                    self.chromosomes[line[0]] = int(line[1])
                else:
                    pass

    def read_TE_file(self, TE_file):
        """
        reads and stores transposon annotation
        """
        with open(TE_file, 'r') as infile:
            for line in infile:
                line = line.rsplit()
                self.transposons[line[4]] = {'chrom':line[0], 'start':line[1], 'stop':line[2]}

    def create_deletion_coords(self, outfile, n):
        """
        randomly select n transposon coordinates and write to outfile
        """
        with open(outfile, 'w+') as dels:
            dels.write("chr\tstart\tstop\tTE\n")
            for x in xrange(n):
                selection = random.choice(self.transposons.keys())
                info = self.transposons[selection]
                dels.write("{ch}\t{sta}\t{stp}\t{TE}\n".format(ch=info['chrom'],
                                                               sta=info['start'],
                                                               stp=info['stop'],
                                                               TE=selection))

    def create_insertion_coords(self, outfile, percCopied, n):
        """
        randomly select n transposon coordinates and randomly generate insertion points
        save to outfile
        """
        n_copied = n*percCopied
        with open(outfile, 'w+') as ins:
            ins.write("chr\tstart\tstop\tdest_chr\tdest_start\tis_copied\n")
            for x in xrange(n):
                selection = random.choice(self.transposons.keys())
                info = self.transposons[selection]
                dest_chrom = random.choice(self.chromosomes.keys())
                dest_start = random.randint(0, self.chromosomes[dest_chrom])
                if n_copied > 0:
                    copied = 'TRUE'
                    n_copied -= 1
                else:
                    copied = 'FALSE'
                ins.write("{ch}\t{sta}\t{stp}\t{dest_ch}\t{dest_sta}\t{copied}\n".format(ch=info['chrom'],
                                                                                         sta=info['start'],
                                                                                         stp=info['stop'],
                                                                                         dest_ch=dest_chrom,
                                                                                         dest_sta=dest_start,
                                                                                         copied=copied))

if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Generates coordinates for synthetic TE insertions / deletions for use in RSVSim')
    parser.add_argument('-g', '--genome', help='genome faidx file', required=True)
    parser.add_argument('-t', '--te', help='TE bedfile', required=True)
    parser.add_argument('-d', '--deletions', help='deletions output file name', required=False, default='deletions.txt')
    parser.add_argument('-i', '--insertions', help='insertions output file name', required=False, default='insertions.txt')
    parser.add_argument('-n', '--number', help='number of transpositions to simulate', required=False, default=100, type=int)
    parser.add_argument('-c', '--copy', help='fraction of transpositions that are copy-paste', required=False, default=0.5, type=float)
    options = parser.parse_args()

    sim = sim_data()
    sim.read_genome_index(options.genome)
    sim.read_TE_file(options.te)
    sim.create_deletion_coords(options.deletions, options.number)
    sim.create_insertion_coords(options.insertions, options.copy, options.number)
