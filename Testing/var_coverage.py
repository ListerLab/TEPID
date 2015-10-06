from subprocess import call
from argparse import ArgumentParser
import os


parser = ArgumentParser(description='Simulate illumina reads of variable coverages from fasta file')
parser.add_argument('-g', '--genome', help='genome fasta file', required=True)
parser.add_argument('-d', '--distance', help='outer distance between the two ends', required=True, type=str)
parser.add_argument('-c', '--coverages', help='comma separated list of coverages to simulate', required=True)
parser.add_argument('-s', '--std', help='standard deviation', required=True, type=str)
parser.add_argument('-p', '--pe', help='paired end reads', default=True, type=bool, required=False)
parser.add_argument('-l', '--length', help='read length', required=True, type=str)
options = parser.parse_args()


def sim_read(options, cov):
    name = cov+'x'
    if options.pe is True:
        multiplier = 2
    else:
        multiplier = 1
    number_read = (int(cov) * 120000000) / (int(options.length) * multiplier)
    if not os.path.exists('./'+name):
        os.makedirs('./'+name)
    call(['wgsim', '-d', options.distance, '-s', options.std, '-N', str(number_read), '-1', options.length, '-2', options.length, options.genome, './{c}/{c}_1.fastq'.format(c=name), './{c}/{c}_2.fastq'.format(c=name)])


coverages = options.coverages.split(',')
for cov in coverages:
    sim_read(options, cov)
