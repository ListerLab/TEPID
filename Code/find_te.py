import pybedtools
import locate
import os


# path to discordant reads bam file
disc_mapped = locate.checkArgs('-d', '--disc')

# path to mapped bam file
all_mapped = locate.checkArgs('-c', '--conc')

# path to split read bam file
split_mapped = locate.checkArgs('-s', '--split')

# path to TE bedfile. Doesn't need to be unzipped
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

# bedtools pairtobed to find TE intersections xor and neither for disc and split reads
te_intersect_disc = disc.pair_to_bed(te, f=0.1, type='xor').moveto('intersect.temp')
no_intersect_disc = disc.pair_to_bed(te, f=0.1, type='neither').moveto('no_intersect.temp')
split_no_intersect = split_bedpe.pair_to_bed(te, f=0.1, type='neither').moveto('split_no_intersect.temp')

# reorder columns so that TE read is in second position
locate.reorder('intersect.temp', 'reorder_intersect.bed')
locate.reorder('no_intersect.temp', 'reorder_no_intersect.bed')
locate.reorder('split_no_intersect.temp', 'reorder_split_no_intersect.bed')

# remove intermediate files
os.remove('intersect.temp')
os.remove('no_intersect.temp')
os.remove('split_no_intersect.temp')

# Now create bedtool objects again, reordered
te_intersect_disc_ordered = pybedtools.BedTool('reorder_intersect.bed').sort()
no_intersect_ordered = pybedtools.BedTool('reorder_no_intersect.bed')
split_no_intersect_ordered = pybedtools.BedTool('reorder_split_no_intersect')
no_intersect_ordered = no_intersect_ordered.cat(split_no_intersect_ordered, postmerge=False).sort()

# merge intersection coordinates where TE name is the same
merged_intersections = locate.merge_TE_coords(te_intersect_disc_ordered, 9)

# filter out where there are reads on different strands
strands = merged_intersections.filter(lambda b: len(b[x]) < 2)  # set x to column containing collapsed strand info
merged_intersections = strands

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
