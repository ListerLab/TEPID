from sys import argv, exit

args = argv[1:]
if '-h' in args or '--help' in args or len(args) == 0:
    print(
            """
            find_te.py

            Created by Tim Stuart

            Usage:
            python find_te.py -n <sample_name>
                              -c <all_mapped_reads_bamfile>
                              -d <discordant_reads_bamfile>
                              -s <split_reads_bamfile>
                              -t <TE_annotation>

            find_te.py must be in the same folder as locate.py
            so that code can be imported

            Outputs TE insertions bedfile and TE deletions bedfile.
            """
            )
    exit()
else:
    pass

import pybedtools
import locate
import merge
import os
from glob import glob
import time
import profile


name = locate.checkArgs('-n', '--name')
all_mapped = locate.checkArgs('-c', '--conc')
# this needs to have unmapped reads removed first due to problem with pybedtools
# samtools view -hbF 0x04 -@ [number_threads] [input] > [output]
print 'Estimating mean insert size'
bam = pybedtools.BedTool(all_mapped)
mn, std = locate.calc_mean(bam[100000:120000])
max_dist = (4*std) + mn

print 'Processing split reads'
split_mapped = locate.checkArgs('-s', '--split')
split = pybedtools.BedTool(split_mapped).bam_to_bed().saveas('split.temp').sort()
breakpoints = split.merge(c='2,3', o='count_distinct,count_distinct')\
.filter(locate.filter_unique_break).each(locate.break_coords).saveas('breakpoints.temp')
locate.convert_split_pairbed('split.temp', 'split_bedpe.temp')
split_bedpe = pybedtools.BedTool('split_bedpe.temp').each(locate.append_origin, word='split').saveas().\
intersect(breakpoints, c=True).each(locate.append_break).saveas().sort()
split_ins = split_bedpe.filter(lambda x: (abs(int(x[1]) - int(x[4])) > 5000) or (x[0] != x[3])).saveas()

# Can't use main bam file because problem when reads don't have their pair,
# due to filtering unmapped reads. This can be fixed when pybedtools problem
# with unmapped reads is fixed
print 'Processing discordant reads'
disc_mapped = locate.checkArgs('-d', '--disc')
# this step takes a long time if there are lots of discordant reads
disc = pybedtools.BedTool(disc_mapped)\
.bam_to_bed(bedpe=True, mate1=True)\
.filter(lambda x: (abs(int(x[1]) - int(x[4])) > max_dist) or (x[0] != x[3])).saveas()\
.each(locate.append_origin, word='disc\tFalse').moveto()
disc_split_dels = split_bedpe.cat(disc, postmerge=False).sort().saveas('disc_split_dels.temp')
disc_split_ins = split_ins.cat(disc, postmerge=False).sort().saveas('disc_split_ins.temp')

print 'Processing TE annotation'
te_bed = locate.checkArgs('-t', '--te')
te = pybedtools.BedTool(te_bed).sort()
te_intersect_ins = disc_split_ins.pair_to_bed(te, f=0.80).saveas('intersect_ins.temp')

print 'Merging TE intersections'
locate.reorder('intersect_ins.temp', 'reorder_intersect.temp')
t0 = time.time()
merge.merge_te_coords('reorder_intersect.temp', 'merged_intersections.temp', 25) 
t1 = time.time()
print (t1-t0)

print 'Finding deletions'
locate.create_deletion_coords(disc_split_dels, 'del_coords.temp')
pybedtools.BedTool('del_coords.temp').intersect(te, wo=True).sort().saveas('deletions.temp')
locate.annotate_deletions('deletions.temp', name, 10, all_mapped, mn)

print 'Finding insertions'
locate.annotate_insertions('merged_intersections.temp', 'insertions.temp')
locate.separate_reads(name)
pybedtools.BedTool('insertions_unsorted.temp').sort().moveto('insertions_{}.bed'.format(name))

temp = glob('./*.temp')
for i in temp:
    os.remove(i)
