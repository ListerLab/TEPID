# note that you can chain pybedtools commands together like unix pipes
# might not need to convert bam to bed as pybedtools can work on bam files
# could also do multiprocessing...

import pybedtools
import locate  # best way to do this?
import os

# need to get input bam file
# assume input is a sorted bam file
# - disc_mapped
# - all_mapped
# - split_mapped
# - te_bed

# create bedtool objects, convert to bed, sort
# need to make object for TE annotation as well
bam = pybedtools.BedTool(all_mapped)
mn, std = locate.calc_mean(bam.head(20000))  # estimate lib size

disc = pybedtools.BedTool(disc_mapped).bam_to_bed(bedpe=True, mate1=True).sort()
split = pybedtools.BedTool(split_mapped).bam_to_bed().sort().saveas('split.temp')
te = pybedtools.BedTool(te_bed).sort()

# remove reads mapped to mito, chlr genomes. May not be needed..

# Convert split bed to bedpe format (already have python script)
# can make a list of strings and create bedtool
# append to same bedtool by using:
# test = pybedtools.BedTool('string', from_string=True).cat(test)

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

Can do this with BedTool.filter() and create a function that will filter out lines
where TE == input name, do for list of unique names. Can go straight from filtered list to
merge, then filter where both strands, then append to final BedTool object.

Can also do parellel_apply
"""

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
