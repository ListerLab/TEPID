import pybedtools
import locate
import os


# path to discordant reads bam file
disc_mapped = locate.checkArgs('-d', '--disc')

# path to mapped bam file
all_mapped = locate.checkArgs('-c', '--conc')

# path to split read bam file
split_mapped = locate.checkArgs('-s', '--split')

# path to TE bedfile
te_bed = locate.checkArgs('-t', '--te')

# create bedtool objects, convert to bed, sort
bam = pybedtools.BedTool(all_mapped)
mn, std = locate.calc_mean(bam.head(20000))  # estimate lib size

disc = pybedtools.BedTool(disc_mapped).bam_to_bed(bedpe=True, mate1=True).sort()
split = pybedtools.BedTool(split_mapped).bam_to_bed().sort().saveas('split.temp')
te = pybedtools.BedTool(te_bed).sort()

# remove reads mapped to mito, chlr genomes. May not be needed..

# Convert split bed to bedpe format (already have python script)
locate.convert_split_pairbed('split.temp', 'split.bedpe')
os.remove('split.temp')
split_bedpe = pybedtools.BedTool('split.bedpe')

# bedtools pairtobed to find TE intersections xor and neither for disc and split reads, sort
te_intersect_disc = disc.pair_to_bed(te, f=0.1, type='xor').sort().moveto('intersect.temp')
no_intersect_disc = disc.pair_to_bed(te, f=0.1, type='neither').sort().moveto('no_intersect.temp')
split_no_intersect = split_bedpe.pair_to_bed(te, f=0.1, type='neither').sort().moveto('split_no_intersect.temp')

# reorder columns so that TE read is in second position
# might be better to create a dictionary and go straight to second part
locate.reorder('intersect.temp', 'reorder_intersect.bed')
locate.reorder('no_intersect.temp', 'reorder_no_intersect.bed')
locate.reorder('split_no_intersect.temp', 'reorder_split_no_intersect.bed')

os.remove('intersect.temp')
os.remove('no_intersect.temp')
os.remove('split_no_intersect.temp')

"""
split intersections file by TE name into different files
define function that reads file, creates bedtool object for each TE name
and does bedtools merge and sort on each bedtool object,
then cats all objects and destroys intermediates,
save final object to file

append to same bedtool by using:
test = pybedtools.BedTool('string', from_string=True).cat(test)

Can do this with BedTool.filter() and create a function that will filter out lines
where TE == input name, do for list of unique names. Can go straight from filtered list to
merge, then filter where both strands, then append to final BedTool object.

Use stream=True to avoid build up of temporary files on disk

Can also do parellel_apply
"""

te_intersect_disc = pybedtools.BedTool('reorder_intersect.bed').sort()


def filter_name(bed, name):
    """
    Filters lines where TE name in bed
    matches input TE name
    """
    for i in bed:
        if bed.name == name:
            print bed
        else:
            pass


# filter out where there are reads on different strands

# cat all merged files

# python annotate_ins.py

# bedtools merge split reads and count unique coords to find breakpoints

# filter out where there is a breakpoint - separate_breakpoints.py

# bedtools intersect for single and double breaks

# python annotate_breakpoints.py
# python separate_reads.py

# sort insertions file and save - use pybedtools moveto (faster)


# Deletions

# python create_deletion_coords.py

# bedtools intersect, python annotate_del.py

# sort deletions file and save
