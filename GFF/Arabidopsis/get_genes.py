from subprocess import call

get_genes = "wget ftp://ftp.arabidopsis.org//Maps/gbrowse_data/TAIR10/TAIR10_GFF3_genes.gff"
get_tes = "wget ftp://ftp.arabidopsis.org//Genes/TAIR9_genome_release/TAIR9_Transposable_Elements.txt"

remove_genes = "rm TAIR10_GFF3_genes.gff filtered_genes.temp TAIR10_genes_temp.bed uniq.temp gene_only_temp.bed"
remove_tes = "rm TAIR9_Transposable_Elements.txt temp_TE.bed"

sort_genes = "sort -k1,1 -nk2,2 filtered_genes.temp > TAIR10_genes_temp.bed"
sort_tes = "sort -k1,1 -nk2,2 temp_TE.bed > TAIR9_TE.bed"

call(get_genes, shell=True)
call(get_tes, shell=True)

with open('TAIR10_GFF3_genes.gff', 'r') as infile, open('filtered_genes.temp', 'w+') as outfile:
    for line in infile:
        line = line.rsplit()
        chrom = line[0]
        feature = line[2]
        if chrom == 'ChrM' or chrom == 'ChrC':
            pass
        elif feature == 'chromosome' or feature == 'CDS':
            pass
        else:
            chrom = chrom.strip('Chr')
            strand = line[6]
            start = line[3]
            stop = line[4]
            info = line[8]
            info = info.split(';')
            AGI = info[0].replace('=', '.').split('.')
            AGI = AGI[1]
            outfile.write("{ch}\t{start}\t{stop}\t{strand}\t{feature}\t{name}\n".format(ch=chrom,
                                                                                        start=start,
                                                                                        stop=stop,
                                                                                        strand=strand,
                                                                                        feature=feature,
                                                                                        name=AGI))

call(sort_genes, shell=True)
call("uniq TAIR10_genes_temp.bed > uniq.temp", shell=True)

with open('uniq.temp', 'r') as infile, open('TAIR10_genes.bed', 'w+') as outfile:
    lines = []
    for i, l in enumerate(infile):
        lines.append(l)
    for x in range(i):
        line = lines[x]
        outfile.write(line)
        line = line.rsplit()
        chrom = line[0]
        feature = line[4]
        stop = line[2]
        strand = line[3]
        AGI = line[5]
        if x+1 <= i and feature == 'exon':
            nextline = lines[x+1]
            nextline = nextline.rsplit()
            nextstart = nextline[1]
            nextfeature = nextline[4]
            nextchrom = nextline[0]
            if nextfeature == 'exon' and nextchrom == chrom and nextstart > stop:
                outfile.write("{ch}\t{start}\t{stop}\t{strand}\t{feature}\t{name}\n".format(ch=chrom,
                                                                                            start=stop,
                                                                                            stop=nextstart,
                                                                                            strand=strand,
                                                                                            feature='intron',
                                                                                            name=AGI))
            else:
                pass
        else:
            pass


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

sort_genes = "sort -k1,1 -nk2,2 filtered_genes.temp > TAIR10_genes_temp.bed"
sort_tes = "sort -k1,1 -nk2,2 temp_TE.bed > TAIR9_TE.bed"

call(sort_tes, shell=True)
call(remove_tes, shell=True)
call("grep 'mRNA' TAIR10_genes.bed > gene_only_temp.bed", shell=True)
call("""awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4,$6}' gene_only_temp.bed > gene_only.bed""", shell=True)
call(remove_genes, shell=True)
