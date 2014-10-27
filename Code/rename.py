import csv
import os
from subprocess import call


with open('SraRunInfo.csv', 'r') as accessions:
    accession_names = {}
    lines = csv.reader(accessions, delimiter=',')
    for line in lines:
        run = line[0]
        name = line[29]
        name = name.replace('(', '_')
        name = name.replace(')', '_')
        accession_names[run] = name

for dirs in os.listdir('.'):
    if dirs in accession_names.keys():
        call("mv {dirs} {acc}".format(dirs=dirs, acc=accession_names[dirs]), shell=True)
        os.chdir(accession_names[dirs])
        for files in os.listdir('.'):
            if files.endswith('_TE_intersections.bed'):
                newname = files.split('_TE')
                acc_name = accession_names[dirs] + '_TE' + newname[1]
            elif files.endswith('.tgz'):
                newname = files.split('.')
                acc_name = accession_names[dirs] + '.' + newname[1]
            elif files.endswith('.fastq'):
                if '_' in files:
                    newname = files.split('_')
                    acc_name = accession_names[dirs] + '_' + newname[1]
                else:
                    newname = files.split('.')
                    acc_name = accession_names[dirs] + '.' + newname[1]
            else:
                newname = files.split('.')
                acc_name = accession_names[dirs] + '.' + newname[1]
            call("mv {files} {new}".format(files=files, new=acc_name), shell=True)
        os.chdir('..')
