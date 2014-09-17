import os
import MySQLdb
import pandas as pd


def get_mc(host, username, password, insertions, context):
    """
    For each insertion, get methylation level in all accessions that have insertion and all that don't have insertion.
    Sum mC for each context in bins of 100 bp for insertion +/- 2kb
    average all bins to get mC state no TE insertion vs mC start TE insertion.
    """
    link = MySQLdb.connect(host, username, password)
    cursor = link.cursor()
    data = {}
    with open(insertions, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            chrom = line[0]
            start = line[1]
            stop = line[2]
            accessions = line[10]
            accessions = accessions.split(',')
            for accession in accessions:
                table = 'mC_calls_{a}'.format(a=accession)
                if start > 2000:
                    upstream = start - 2000
                    if end < chr_end:
                        downstream = stop + 2000
                        for x in range(40):
                            bins = 100
                            bins_start = (upstream + (bins * (x-1)))
                            bins_end = (upstream + (bins * x))
                            query = """select sum(mc), sum(h) from population_epigenetics.{table}
                                        where class = '{context}' assembly = {ch} and (position between {start} and {end})""".format(table=table,
                                                                                                                                     ch=chrom,
                                                                                                                                     context=context,
                                                                                                                                     start=bins_start,
                                                                                                                                     end=bins_end)
                            # upstream
                            if x < 20:
                                addDataC(query, data, x, cursor)
                            # Downstream
                            elif x >= 20:
                                addDataC(query, data, x, cursor)
                            else:
                                pass
                    else:
                        pass
                else:
                    pass
        error = {}
        for key, value in data.iteritems():
            if sum(value) > 0:
                mean = sum(value)/len(value)
            else:
                mean = 0.0
            data[key] = round(mean, 5)
        cursor.close()
        link.close()
        return data


def addDataC(q, d, x, cursor):
    cursor.execute(q)
    results = cursor.fetchall()
    for row in results:
        mc = row[0]
        h = row[1]
    if h == 0 or h is None:  # no coverage
        pass
    elif h > 0:
        level = mc / h
        level = float(level)
        try:
            d[x]
        except KeyError:
            d[x] = [level]
        else:
            d[x].append(level)


def saveData(filename, *args):
    """
    Creates pandas dataframe and appends each dictionary in args as a column.
    Saves dataframe to csv file.
    """
    l = len(args)
    df = pd.DataFrame.from_dict(args[0], orient='index', dtype=None)
    df.columns = [args[1]]
    for x in range(2, l-1, 2):
        if args[x] is not None:
            df[args[x+1]] = pd.DataFrame.from_dict(args[x], orient='index', dtype=None)
        else:
            pass
    df.to_csv(filename, sep='\t')


data_cg = get_mc(host, username, password, 'sorted_insertions.bed', 'CG')
data_chg = get_mc(host, username, password, 'sorted_insertions.bed', 'CHG')
data_chh = get_mc(host, username, password, 'sorted_insertions.bed', 'CHH')

saveData('mc_insertions.tsv',
         data_cg, 'cg',
         data_chg, 'chg',
         data_chh, 'chh')
