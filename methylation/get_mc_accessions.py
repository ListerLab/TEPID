import os

with open('c_dmrs.tsv', 'r') as dmr_file:
    dmrs = {}
    lines = []
    for i, l in enumerate(dmr_file):
        lines.append(l)
    for x in range(i):
        line = lines[x]
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
            if line[0] == 'chrom':
                pass  # header
            else:
                methylation_call = line[6]
                position = line[1]
                chromosome = line[0]
                if value['chrom'] > chromosome:  # might make it a bit faster
                    pass
                elif value['start'] < position < value['stop'] and chromosome == value['chrom']:
                    total_mc += int(methylation_call)
                elif value['stop'] < position:
                    break
                else:
                    pass
        try:
            levels[key]
        except KeyError:
            levels[key] = [total_mc]
        else:
            levels[key].append(total_mc)

files = [f for f in os.listdir('.') if os.path.isfile(f)]
levels = {}
header = []

for f in files:
    if f.startswith('GSM'):
        with open(f, 'r') as infile:
            accession_name = f.replace('calls_', '.').split('.')
            accession_name = accession_name[1]
            header.append(accession_name)
            add_data(accession_name, infile)

with open('dmrs_acc.tsv', 'w+') as outfile:
    outfile.write('chr\tstart\tstop\t{a}\n'.format(a='\t'.join(header)))
    for key, value in dmrs.items():
        outfile.write("{chrom}\t{start}\t{stop}\t{a}\n".format(chrom=value["chrom"],
                                                               start=value["start"],
                                                               stop=value["stop"],
                                                               a='\t'.join(map(str, levels[key]))))
