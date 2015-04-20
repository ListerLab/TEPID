from subprocess import call


def sim_read(cov, read_length, d, s, name, genome, pair=True):
    d = str(d)
    s = str(s)
    if pair is True:
        multiplier = 2
    else:
        multiplier = 1
    number_read = (cov * 120000000) / (read_length * multiplier)
    call(['wgsim', '-d', d, '-s', s, '-N', str(number_read), '-1', str(read_length), '-2', str(read_length), genome, '{}_1.fastq'.format(name), '{}_2.fastq'.format(name)])


# 5x
sim_read(5, 100, 250, 50, '5x', 'genome_rearranged.fasta', pair=True)

# 10x
sim_read(10, 100, 250, 50, '10x', 'genome_rearranged.fasta', pair=True)

# 20x
sim_read(20, 100, 250, 50, '20x', 'genome_rearranged.fasta', pair=True)

# 30x
sim_read(30, 100, 250, 50, '30x', 'genome_rearranged.fasta', pair=True)

# 40x
sim_read(40, 100, 250, 50, '40x', 'genome_rearranged.fasta', pair=True)

# 50x
sim_read(50, 100, 250, 50, '50x', 'genome_rearranged.fasta', pair=True)

# 60x
sim_read(60, 100, 250, 50, '60x', 'genome_rearranged.fasta', pair=True)

# 70x
sim_read(70, 100, 250, 50, '70x', 'genome_rearranged.fasta', pair=True)