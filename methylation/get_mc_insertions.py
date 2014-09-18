import os
import MySQLdb
import pandas as pd
from sys import argv


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


def get_mc(host, username, password, insertions, context):
    """
    For each insertion, get methylation level in all accessions that have insertion and all that don't have insertion.
    Sum mC for each context in bins of 100 bp for insertion +/- 2kb
    average all bins to get mC state no TE insertion vs mC start TE insertion.
    """
    link = MySQLdb.connect(host, username, password)
    cursor = link.cursor()
    data = {}
    cursor.execute("use population_epigenetics;")
    cursor.execute("show tables")
    tables_result = cursor.fetchall()
    names = []
    for row in tables_result:
        name = row[0]
        names.append(name)
    with open(insertions, 'r') as infile:
        for line in infile:
            line = line.rsplit()
            chrom = int(line[0])
            start = int(line[1])
            stop = int(line[2])
            insertion_point = int((stop+start) / 2)
            accessions = line[10]
            accessions = accessions.split(',')
            for accession in accessions:
                accession = accession.replace('-', '_')
                print accession
                table = 'mC_calls_{a}'.format(a=accession)
                if table in names:
                    # get chromosome end
                    chr_query = "select position from {tb} where assembly = 1 order by position desc limit 1".format(tb=table)
                    cursor.execute(chr_query)
                    chr_results = cursor.fetchall()
                    for row in chr_results:
                        chr_end = int(row[0])
                    if start > 2000:
                        upstream = insertion_point - 2000
                        for x in range(40):
                            bins = 100
                            bins_start = (upstream + (bins * (x-1)))
                            bins_end = (upstream + (bins * x))
                            query = """select sum(mc), sum(h) from {table} where class = '{context}'
                                       and assembly = {ch} and (position between {start} and {end})""".format(table=table,
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

host = checkArgs('h', 'host')
username = checkArgs('u', 'username')
password = checkArgs('p', 'password')
infile = checkArgs('f', 'file')

data_cg = get_mc(host, username, password, infile, 'CG')
data_chg = get_mc(host, username, password, infile, 'CHG')
data_chh = get_mc(host, username, password, infile, 'CHH')

saveData('mc_insertions.tsv',
         data_cg, 'cg',
         data_chg, 'chg',
         data_chh, 'chh')
