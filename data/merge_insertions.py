import os


def overlap(start1, stop1, start2, stop2):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    for y in xrange(start2, stop2):
        if start1 <= y <= stop1:
            return True
        else:
            pass


def create_master_dict(master, accession_name):
    with open(master, 'r') as masterfile:
        master_insertions = {}
        for line in masterfile:
            line = line.rsplit()
            if line[0] == 'ins_chrom':
                pass
            else:
                master_insertions[line[5]] = {'ins_chrom': line[0], 'ins_start': line[1], 'ins_end': line[2],
                                              'ins_strand': line[3], 'agi': line[4], 'ref_chr': line[6],
                                              'ref_start': line[7], 'ref_end': line[8], "ref_strand": line[9], 'accessions': [accession_name]}


def merge_insertions(master_dict, ins_file, accession_name):
    with open(ins_file, 'r') as insertions:
        for line in insertions:
            line = line.rsplit()
            ins_chrom = line[0]
            ins_start = line[1]
            ins_end = line[2]
            agi = line[4]
            for key, value in master_dict.items():
                if value['ins_chrom'] == ins_chrom and value['agi'] == agi:
                    if overlap(value['ins_start'], value['ins_end'], ins_start, ins_end) is True:
                        # need to do something with unique id
                        value['accessions'].append(accession_name)
                        # refine insertion coordinates
                        if value['end'] > ins_start > value['ins_start']:
                            value['ins_start'] = ins_start
                        elif value['start'] < ins_end < value['ins_end']:
                            value['ins_end'] = ins_end
                        else:
                            pass
                    else:
                        pass
                else:
                    pass


for dirs in os.listdir('.'):
    if os.path.isdir(dirs) is True:
        os.chdir(dirs)
        if os.path.isfile('insertions_{d}.bed'.format(d=dirs)) is True:
            print dirs
            try:
                master_insertions
            except NameError:
                create_master_dict('insertions_{d}.bed'.format(d=dirs), dirs)
            else:
                merge_insertions(master_insertions, 'insertions_{d}.bed'.format(d=dirs))
            os.chdir('..')
        else:
            os.chdir('..')
    else:
        pass


# save master insertions dict to csv file
with open('insertions.bed', 'w+') as outfile:
    outfile.write("ins_chrom\tins_start\tins_end\tins_strand\tagi\tref_chr\tref_start\tref_end\tref_strand\taccessions\n")
    for key, value in master_insertions.items():
        outfile.write("""{ins_chrom}\t{ins_start}\t{ins_end}\t{ins_strand}\t{agi}\t
                         {ref_chrom}\t{ref_start}\t{ref_end}\t{ref_strand}\t{id}\t{accessions}\n""".format(ins_chrom=value['ins_chrom'],
                                                                                                           ins_start=value['ins_start'],
                                                                                                           ins_end=value['ins_end'],
                                                                                                           ins_strand=value['ins_strand'],
                                                                                                           agi=value['agi'],
                                                                                                           ref_chrom=value['ref_chom'],
                                                                                                           ref_end=value['ref_end'],
                                                                                                           ref_strand=value['ref_strand'],
                                                                                                           id=value['id'],
                                                                                                           accessions=','.join(map(str, value['accessions']))))
