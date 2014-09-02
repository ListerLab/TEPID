from subprocess import call


with open('SraAccList.txt', 'r') as accessions:
    for row in accessions:
        acc = row.strip('\n')
        print "Downloading {a}".format(a=acc)
        call(["wget",
              "-r",
              "--no-parent",
              "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR492/{acc}/*".format(acc=acc)])
