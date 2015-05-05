#! /usr/local/bin/python

import random
from sys import argv, exit

args = argv[1:]
if '-h' in args or '--help' in args or len(args) == 0:
    print(
            """
            generate_sim_coords.py

            Created by Tim Stuart

            Generates coordinates for synthetic TE insertions / deletions
            for use in RSVSim package

            Usage:
            python generate_sim_coords.py -g <genome.fa.fai>
                                          -t <te_bedfile>
                                          -d <deletions_filename>
                                          -i <insertions_filename>
            """
            )
    exit()
else:
    pass

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
                else:
                    self.chromosomes[line[0]] = int(line[1])

    def read_TE_file(self, TE_file):
        """
        reads and stores transposon annotation
        """
        with open(TE_file, 'r') as infile:
            for line in infile:
                line = line.rsplit()
                self.transposons[line[4]] = {'chrom':line[0], 'start':line[1], 'stop':line[2]}

    def create_deletion_coords(self, outfile, n=100):
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

    def create_insertion_coords(self, outfile, percCopied=0.5, n=100):
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
    genome = checkArgs('-g', '--genome')
    te = checkArgs('-t', '--te')
    del_file = checkArgs('-d', '--deletions')
    ins_file = checkArgs('-i', '--insertions')

    sim = sim_data()
    sim.read_genome_index(genome)
    sim.read_TE_file(te)
    sim.create_deletion_coords(del_file)
    sim.create_insertion_coords(ins_file)
