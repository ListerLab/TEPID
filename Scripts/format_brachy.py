# -*- coding: utf-8 -*-


def reformat(fname, out):
    """
    Reformat Brachypodium gff3 TE annotation file to bedfile
    """
    conversion = {'RLC': 'LTR/Copia', 'RLG': 'LTR/Gypsy', 'RLB': 'LTR/Bel–Pao', 'RLR': 'LTR/Retrovirus',
                  'RLE': 'LTR/ERV', 'RLX': 'LTR/unknown', 'RYD': 'DIRS/DIRS', 'RYN': 'DIRS/Ngaro', 'RYV': 'DIRS/VIPER',
                  'RYX': 'DIRS/unknown', 'RPP': 'PLE/Penelope', 'RPX': 'PLE/unknown', 'RIR': 'LINE/R2', 'RIT': 'LINE/RTE', 'RIJ': 'LINE/Jockey',
                  'RIL': 'LINE/L1', 'RII': 'LINE/I', 'RIX': 'LINE/unknown', 'RST': 'SINE/tRNA', 'RSL': 'SINE/7SL', 'RSS': 'SINE/5S', 'RSX': 'SINE/unknown',
                  'DTT': 'TIR/Tc1-Mariner', 'DTA': 'TIR/hAT', 'DTM': 'TIR/Mutator', 'DTE': 'TIR/Merlin',
                  'DTR': 'TIR/Transib', 'DTP': 'TIR/P', 'DTB': 'TIR/PiggyBac', 'DTH': 'TIR/PIF-Harbinger',
                  'DTC': 'TIR/CACTA', 'DTX': 'TIR/unknown', 'DYC': 'Crypton/Crypton', 'DYX': 'Crypton/unknown', 'DHH': 'Helitron/Helitron',
                  'DHX': 'Helitron/unknown', 'DMM': 'Maverick/Maverick', 'DMX': 'Maverick/unknown'}
    with open(fname, 'r') as infile, open(out, 'w+') as outfile:
        for _ in xrange(18):  # skips header lines
            next(infile)
        for line in infile:
            line = line.rsplit()
            raw_chrom = line[0]
            if 'scaffold' in raw_chrom or '##' in raw_chrom:
                pass
            else:
                chrom = raw_chrom.strip('Bd')
                start = line[3]
                stop = line[4]
                strand = line[6]
                info = line[8]
                te_type = line[2]
                info = info.split(';')
                if te_type == 'transposable_element':
                    temp = {}
                    for item in info:
                        splits = item.split('=')
                        temp[splits[0]] = splits[1]
                    ident = temp['ID']
                    name = temp['Name']
                    name = name[4:]
                    name = name.split('_Bd{ch}'.format(ch=chrom))
                    name = name[0]
                    te_class_key = temp['class']
                    te_class = conversion[te_class_key]
                    te_id = 'Bd{ch}TE{id}'.format(ch=chrom, id=ident)
                    outfile.write('{ch}\t{sta}\t{sto}\t{str}\t{id}\t{nm}\t{cls}\n'.format(ch=chrom, sta=start, sto=stop, str=strand, id=te_id, nm=name, cls=te_class))
                else:
                    pass

if __name__ == '__main__':
    reformat('Brachy_TEs_v2.2_16-07-2009.gff3', 'Brachy_TE_v2.2.bed')
