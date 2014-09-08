import os

with open('c_dmrs.tsv', 'r') as dmr_file:
    dmrs = {}
    num_lines = sum(1 for line in dmr_file)
    for x in range(num_lines):
        line = dmr_file[x]
        line = line.rsplit()
        chrom = line[0].strip('chr')
        start = line[1]
        stop = line[2]
        dmrs[x] = {'chrom': chrom, 'start': start, 'stop': stop}


def add_data(accession, acc_file):
    for key, value in dmrs.items():
        total_mc = 0
        for line in acc_file:
            line = line.rsplit()
            methylated_bases =  # find position of this value
            position = 
            chromosome = 
            if value['start'] < postion > value['stop'] and chromosome == value['chrom']:
                total_mc += int(methylated_bases)
            else:
                pass
        levels.append(total_mc)
        header.append(accession)  # keeps order correct


files = [f for f in os.listdir('.') if os.path.isfile(f)]
levels = []
header = []

for f in files:
    with open(f, 'r') as infile:
        accession_name = f.replace('calls_', '.')
        accession_name = f.split('.')
        accession_name = accession_name[1]
        add_data(accession_name, infile)


with open('dmrs_acc.tsv', 'w+') as outfile:
    outfile.write('chr\tstart\tstop\t{a}\n'.format(a='\t'.join(header)))
    for key, value in dmrs.items():
        outfile.write("{chrom}\t{start}\t{stop}\t{a}\n".format(chrom=value["chrom"],
                                                               start=value["start"],
                                                               stop=value["stop"],
                                                               a='\t'.join(levels)))
