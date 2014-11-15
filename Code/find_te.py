# note that you can chain pybedtools commands together like unix pipes
# might not need to convert bam to bed as pybedtools can work on bam files
# could also do multiprocessing...

import pybedtools
import locate  # best way to do this?

# need to get input bam file
# assume input is a sorted bam file
# - disc_mapped
# - split_mapped
# - te_bed

# create bedtool objects, convert to bed, sort
# need to make object for TE annotation as well
disc_bam = pybedtools.BedTool(disc_mapped)
mn, std = locate.calc_mean(disc_bam.head(20000))  # estimate lib size

disc = disc_bam.bam_to_bed(bedpe=True, mate1=True).sort()
split = pybedtools.BedTool(split_mapped).bam_to_bed().sort()
te = pybedtools.BedTool(te_bed).sort()

# remove reads mapped to mito, chlr genomes. May not be needed..

# Convert split bed to bedpe format (already have python script)
split_bedpe = locate.convert_split_pairbed(split_bed)  # script needs to be edited

# bedtools pairtobed to find TE intersections xor and neither for disc and split reads, sort
te_intersect_disc = disc.pair_to_bed(te, f=0.1, type='xor')
no_intersect_disc = disc.pair_to_bed(te, f=0.1, type='neither')

# reorder columns so that TE read is in second position (reorder.py does this, need to update code)
locate.reorder(te_intersect_disc)
locate.reorder(no_intersect_disc)

# sort by TE name and split into file / objects for each TE

# bedtools sort and merge


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
