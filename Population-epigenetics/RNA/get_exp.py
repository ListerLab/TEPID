import numpy as np


def get_exp(genes, fpkms):
    """
    Create dictionary with expression data for all genes for accession.
    Sum expression of genes with TE insertions.
    Genes parameter is set of genes to gather expression data for.
    FPKMs parameter is file containing expression data.
    """
    with open(genes, 'r') as accession, open(fpkms, 'r') as expression:
        names = []
        exp = {}
        total = []
        for line in accession:
            line = line.rsplit()
            name = line[8]
            if name == '.':
                pass
            else:
                names.append(name)
        for line in expression:
            line = line.rsplit()
            if line[0] == 'tracking_id':
                pass
            else:
                agi = line[0]
                fpkm = line[2]
                if 'e' in fpkm:
                    fpkm = fpkm.split('e')
                    multiplication = fpkm[1]
                    base = float(fpkm[0])
                    if '+' in multiplication:
                        multiplication = multiplication.strip('+')
                        fpkm = base * 10**int(multiplication)
                        exp[agi] = fpkm
                    elif '-' in multiplication:
                        multiplication = multiplication.strip('-')
                        fpkm = base * 10**(-int(multiplication))
                        exp[agi] = fpkm
                    else:
                        raise Exception("error: invalid fpkm values")
                else:
                    fpkm = float(fpkm)
                    exp[agi] = fpkm
        for name in names:
            total.append(exp[name])
    add = sum(total)
    number = len(total)
    mean = add / number
    std = np.std(total)
    rtn = np.sqrt(number)
    sem = std / rtn
    return mean, sem

exp_list = []

exp_list.append(get_exp('upstream_insertions.bed', 'Aa_0_expression.tsv'))
exp_list.append(get_exp('upstream_insertions.bed', 'Col_0_expression.tsv'))

exp_list.append(get_exp('gene_overlap.bed', 'Aa_0_expression.tsv'))
exp_list.append(get_exp('gene_overlap.bed', 'Col_0_expression.tsv'))

exp_list.append(get_exp('downstream_insertions.bed', 'Aa_0_expression.tsv'))
exp_list.append(get_exp('downstream_insertions.bed', 'Col_0_expression.tsv'))

newlist = []
for item in exp_list:
    newlist.append(str(item[0]))
    newlist.append(str(item[1]))

with open('insertion_expression.tsv', 'w+') as outfile:
    outfile.write('up_aa0\tup_aa0_sem\tup_col0\tup_col0_sem\tin_aa0\tin_aa0_sem\tin_col0\tin_col0_sem\tdn_aa0\tdn_aa0_sem\tdn_col0\tdn_col0_sem\n')
    outfile.write('{a}\n'.format(a='\t'.join(newlist)))
