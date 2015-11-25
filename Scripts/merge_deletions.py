#! /usr/bin/env python


def create_master_dict(sample, fname):
    with open(fname, 'r') as masterfile:
        x = 0
        master_dict = {}
        for line in masterfile:
            field = line.rsplit()
            if line[0] == 'ins_chr':
                pass
            else:
                coords = '\t'.join(field[:5])
                master_dict[x] = {'coords': coords, 'accessions': [sample]}
                x += 1
        return master_dict


def merge_deletions(master, fname, sample):
    with open(fname, 'r') as infile:
        for line in infile:
            field = line.rsplit()
            coords = '\t'.join(field[:5])
            i = len(master)-1
            x = 0
            while x <= i:
                if master[x]['coords'] == coords:
                    master[x]['accessions'].append(sample)
                    break
                elif x == i:
                    master[x+1] = {'coords': coords, 'accessions': [sample]}
                    break
                else:
                    x += 1


def save_deletions(master, outf):
    with open(outf, 'w+') as outfile:
        for key, value in master.iteritems():
            accessions = set(value['accessions'])
            outfile.write('{c}\t{a}\n'.format(c=value['coords'], a=','.join(accessions)))

if __name__ == "__main__":
    import os
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Merge TE deletions calls')
    parser.add_argument('-f', '--filename', help='filename prefix for merge files', required=True)
    options = parser.parse_args()

    for dirs in os.listdir('.'):
        if os.path.isdir(dirs) is True:
            os.chdir(dirs)
            if os.path.isfile(options.filename+'_{d}.bed'.format(d=dirs)) is True:
                print "Processing {d}".format(d=dirs)
                try:
                    master_dictionary
                except NameError:
                    master_dictionary = create_master_dict(dirs, options.filename+'_{d}.bed'.format(d=dirs))
                else:
                    merge_deletions(master_dictionary, options.filename+'_{d}.bed'.format(d=dirs), dirs)
                os.chdir('..')
            else:
                os.chdir('..')
        else:
            pass
    save_deletions(master_dictionary, options.filename+'.bed')
