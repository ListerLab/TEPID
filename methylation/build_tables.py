# usage:
# python build_tables.py h <host> u <username> p <password> d <database>

from sys import argv
import MySQLdb


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

"""
MySQL commands:
create database population_epigenetics;
use population_epigenetics;
create table mC_calls_Aa_0 (assembly int(1), position int(9), strand char(1), class char(3), mc int(3), h int(3), uc int(3));
insert into mC_calls_Aa_0 (assembly, position, strand, class, mc, h, uc) values (1, 10, '+', 'CHH', 4, 5, 1);
"""


def add_table(host, username, password, database, fname):
    """Make new table and add data"""
    with open(fname, 'r') as infile:
        link = MySQLdb.connect(host, username, password)
        cursor = link.cursor()
        accession_name = fname.replace('calls_', '.').split('.')
        accession_name = accession_name[1]
        cursor.execute("use {db};".format(db=database))
        build = """create table mC_calls_{name} (assembly int(1), position int(9), strand char(1), class char(3), mc int(3), h int(3), uc int(3));""".format(name=accession_name)
        cursor.execute(build)
        for line in infile:
            line = line.rsplit()
            chrom = line[0]
            if chrom == 'chrom':
                pass  # header
            else:
                chrom = int(chrom)
                pos = int(line[1])
                strand = line[2]
                mc_class = line[3]
                mc_class = mc_class[0] + mc_class[1:].replace('A', 'H').replace('T', 'H').replace('C', 'H')
                mc = int(line[4])
                h = int(line[5])
                uc = h - mc
                add_data = """insert into mC_calls_{name} (assembly, position, strand, class, mc, h, uc) values ({chrom}, {pos}, {strand}, {class}, {mc}, {h}, {uc});""".format(name=name,
                                                                                                                                                                                chrom=chrom,
                                                                                                                                                                                pos=pos,
                                                                                                                                                                                strand=strand,
                                                                                                                                                                                class=mc_class,
                                                                                                                                                                                mc=mc,
                                                                                                                                                                                h=h,
                                                                                                                                                                                uc=uc)
                cursor.execute(add_data)
        cursor.close()
        link.close()


host = checkArgs(h, host)
username = checkArgs(u, user)
password = checkArgs(p, password)
database = checkArgs(d, database)

for files in os.listdir('.'):
    if files.endswith(".tsv"):
        add_table(host, username, password, database, files)
