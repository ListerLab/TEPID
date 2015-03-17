from subprocess import call

get_tes = "wget ftp://ftp.arabidopsis.org//Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt"
remove_tes = "rm TAIR9_Transposable_Elements.txt temp_TE.bed"
sort_tes = "sort -k1,1 -nk2,2 temp_TE.bed > TAIR9_TE.bed"
call(get_tes, shell=True)


with open('TAIR9_Transposable_Elements.txt', 'r') as infile, open('temp_TE.bed', 'w+') as outfile:
    for line in infile:
        line = line.rsplit()
        name = line[0]
        if name == 'Transposon_Name':
            pass  # skips header
        else:
            strand = line[1]
            start = line[2]
            stop = line[3]
            family = line[4]
            superfamily = line[5]
            chrom = name.split('T')
            chrom = chrom[1]
            if strand == 'false':
                strand = '-'
            elif strand == 'true':
                strand = '+'
            else:
                raise Exception('Check strand information')
            outfile.write("{ch}\t{start}\t{stop}\t{strand}\t{name}\t{family}\t{superfamily}\n".format(ch=chrom,
                                                                                                      start=start,
                                                                                                      stop=stop,
                                                                                                      strand=strand,
                                                                                                      name=name,
                                                                                                      family=family,
                                                                                                      superfamily=superfamily))

sort_tes = "sort -k1,1 -nk2,2 temp_TE.bed > TAIR9_TE.bed"

call(sort_tes, shell=True)
call(remove_tes, shell=True)
