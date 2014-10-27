
def separate_reads(acc):
    with open('insertions_{a}_temp.bed'.format(a=acc), 'r') as infile, open('insertions_{a}_unsorted.bed'.format(a=acc), 'w+') as outfile, open('id_{a}.fa'.format(a=acc), 'w+') as id_file:
        x = 0
        for line in infile:
            line = line.rsplit()
            data = line[:9]
            reads = line[9]
            mates = line[10]
            ident = acc + '_' + str(x)
            outfile.write('{data}\t{id}\n'.format(data='\t'.join(data), id=ident))
            id_file.write('>{id}\t{reads}\t{mates}\n'.format(id=ident, reads=reads, mates=mates))
            x += 1

if __name__ == '__main__':
    import sys
    separate_reads(sys.argv[1])
