from subprocess import call
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

path = checkArgs('-p', '--path')

with open(path, 'r') as accessions:
    for row in accessions:
        acc = row.strip('\n')
        print "Downloading {a}".format(a=acc)
        call("wget -r --no-parent ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR492/{acc}/*".format(acc=acc), shell=True)
