import pybedtools
import locate
import os
from glob import glob


# sample name
name = locate.checkArgs('-n', '--name')

# mapped bam file
# all_mapped = locate.checkArgs('-c', '--conc')
# bam = pybedtools.BedTool(all_mapped)
# mn, std = locate.calc_mean(bam.head(20000))  # estimate lib size. This gives an error (I think because of the header)

# split read data
print 'Processing split reads'
split_mapped = locate.checkArgs('-s', '--split')
split = pybedtools.BedTool(split_mapped).bam_to_bed().sort().saveas('split.bed')
locate.convert_split_pairbed('split.bed', 'split.bedpe')
split_bedpe = pybedtools.BedTool('split.bedpe')
merged_split = split.merge(c='2,3', o='count_distinct,count_distinct')

# need to filter out mito, chlr genes before this step
locate.filter_split(merged_split)  # error loading bedtool file midway through processing
filtered_split = pybedtools.BedTool('filtered_split.temp')

# discordant reads
print 'Processing discordant reads'
disc_mapped = locate.checkArgs('-d', '--disc')
disc = pybedtools.BedTool(disc_mapped).bam_to_bed(bedpe=True, mate1=True).sort()

# TE bedfile
print 'Processing TE annotation'
te_bed = locate.checkArgs('-t', '--te')
te = pybedtools.BedTool(te_bed).sort()

# bedtools pairtobed to find TE intersections xor and neither for disc and split reads
print 'Intersecting TE coordinates with reads'
te_intersect_disc = disc.pair_to_bed(te, f=0.1, type='xor').moveto('intersect.temp')
no_intersect_disc = disc.pair_to_bed(te, f=0.1, type='neither')
no_intersect_split = split_bedpe.pair_to_bed(te, f=0.1, type='neither')
no_intersect = no_intersect_disc.cat(no_intersect_split, postmerge=False).sort()

# reorder columns so that TE read is in second position
locate.reorder('intersect.temp', 'reorder_intersect.bed')
# should filter for discordant reads here to reduce file sizes in later steps

# Now create bedtool objects again, reordered
te_intersect_disc_ordered = pybedtools.BedTool('reorder_intersect.bed').sort()

# merge intersection coordinates where TE name is the same
print 'Merging TE intersections'
merged_intersections = locate.merge_TE_coords(te_intersect_disc_ordered, 9)  # still very slow (>1 hour). Multiprocessing? Also reduce file size before this step

print "Finished merging intersections. Filtering strands"
# filter out where there are reads on different strands
strands = merged_intersections.filter(lambda b: len(b[3]) < 2).saveas('collapsed_intersections.temp')

# Find where there are breakpoints at insertion site
print 'Annotating insertions'
locate.annotate_insertions('collapsed_intersections.temp', 'insertions.temp', name, mn, std)  # mn, sdt not working now due to problem with bam file
pybedtools.BedTool('insertions.temp').intersect(filtered_split, c=True).moveto('insertions_split_reads.temp')  # problem where chromosomes in one file have 'chr', the other doesn't (reports no intersections)
locate.splitfile('insertions_split_reads.temp')  # separate breakpoints. make single and double breakpoint files
pybedtools.BedTool('single_break.temp').intersect(filtered_split, wo=True).moveto('single_break_intersect.temp')
pybedtools.BedTool('double_break.temp').intersect(filtered_split, wo=True).moveto('double_break_intersect.temp')
locate.annotate_single_breakpoint()
locate.annotate_double_breakpoint()
locate.separate_reads(name)
pybedtools.BedTool('insertions_unsorted.temp').sort().moveto('insertions_{a}.bed'.format(a=name))


# Deletions

# python create_deletion_coords.py

# bedtools intersect, python annotate_del.py

# sort deletions file and save

# remove temp files
# temp = glob('./*.temp')
# for i in temp:
#     os.remove(i)

