library(RSVSim)

# load genomes
tair10 = readDNAStringSet("tair10.fa")
brachy = readDNAStringSet("Bdistachyon_283_assembly_v2.0.fa")

# load data and read into correct format
del_tair_r1 <- read.table("tair_del_r1.bed", header=TRUE)
del_tair_r2 <- read.table("tair_del_r2.bed", header=TRUE)
del_tair_r3 <- read.table("tair_del_r3.bed", header=TRUE)
del_tair_r4 <- read.table("tair_del_r4.bed", header=TRUE)
del_tair_r5 <- read.table("tair_del_r5.bed", header=TRUE)

ins_tair_r1 <- read.table("tair_ins_r1.bed", header=TRUE)
ins_tair_r2 <- read.table("tair_ins_r2.bed", header=TRUE)
ins_tair_r3 <- read.table("tair_ins_r3.bed", header=TRUE)
ins_tair_r4 <- read.table("tair_ins_r4.bed", header=TRUE)
ins_tair_r5 <- read.table("tair_ins_r5.bed", header=TRUE)

del_brachy_r1 <- read.table("brachy_del_r1.bed", header=TRUE)
del_brachy_r2 <- read.table("brachy_del_r2.bed", header=TRUE)
del_brachy_r3 <- read.table("brachy_del_r3.bed", header=TRUE)
del_brachy_r4 <- read.table("brachy_del_r4.bed", header=TRUE)
del_brachy_r5 <- read.table("brachy_del_r5.bed", header=TRUE)

ins_brachy_r1 <- read.table("brachy_ins_r1.bed", header=TRUE)
ins_brachy_r2 <- read.table("brachy_ins_r2.bed", header=TRUE)
ins_brachy_r3 <- read.table("brachy_ins_r3.bed", header=TRUE)
ins_brachy_r4 <- read.table("brachy_ins_r4.bed", header=TRUE)
ins_brachy_r5 <- read.table("brachy_ins_r5.bed", header=TRUE)

ath_deletions_r1 <- with(del_tair_r1, GRanges(chr, IRanges(start, stop)))
ath_deletions_r2 <- with(del_tair_r2, GRanges(chr, IRanges(start, stop)))
ath_deletions_r3 <- with(del_tair_r3, GRanges(chr, IRanges(start, stop)))
ath_deletions_r4 <- with(del_tair_r4, GRanges(chr, IRanges(start, stop)))
ath_deletions_r5 <- with(del_tair_r5, GRanges(chr, IRanges(start, stop)))

brachy_deletions_r1 <- with(del_brachy_r1, GRanges(chr, IRanges(start, stop)))
brachy_deletions_r2 <- with(del_brachy_r2, GRanges(chr, IRanges(start, stop)))
brachy_deletions_r3 <- with(del_brachy_r3, GRanges(chr, IRanges(start, stop)))
brachy_deletions_r4 <- with(del_brachy_r4, GRanges(chr, IRanges(start, stop)))
brachy_deletions_r5 <- with(del_brachy_r5, GRanges(chr, IRanges(start, stop)))

ath_insertions_r1 <- with(ins_tair_r1, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r2 <- with(ins_tair_r2, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r3 <- with(ins_tair_r3, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r4 <- with(ins_tair_r4, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
ath_insertions_r5 <- with(ins_tair_r5, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))

brachy_insertions_r1 <- with(ins_brachy_r1, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
brachy_insertions_r2 <- with(ins_brachy_r2, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
brachy_insertions_r3 <- with(ins_brachy_r3, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
brachy_insertions_r4 <- with(ins_brachy_r4, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))
brachy_insertions_r5 <- with(ins_brachy_r5, GRanges(chr, IRanges(start, stop), chrB=dest_chr, startB=dest_start, copied=is_copied))

# simulate variation
simulateSV(output="./ath/r1",
           genome=tair10,
           regionsDels=ath_deletions_r1,
           regionsIns=ath_insertions_r1,
           random=FALSE)

simulateSV(output="./ath/r2",
           genome=tair10,
           regionsDels=ath_deletions_r2,
           regionsIns=ath_insertions_r2,
           random=FALSE)

simulateSV(output="./ath/r3",
           genome=tair10,
           regionsDels=ath_deletions_r3,
           regionsIns=ath_insertions_r3,
           random=FALSE)

simulateSV(output="./ath/r4",
           genome=tair10,
           regionsDels=ath_deletions_r4,
           regionsIns=ath_insertions_r4,
           random=FALSE)

simulateSV(output="./ath/r5",
           genome=tair10,
           regionsDels=ath_deletions_r5,
           regionsIns=ath_insertions_r5,
           random=FALSE)


simulateSV(output="./brachy/r1",
           genome=brachy,
           regionsDels=brachy_deletions_r1,
           regionsIns=brachy_insertions_r1,
           random=FALSE)

simulateSV(output="./brachy/r2",
           genome=brachy,
           regionsDels=brachy_deletions_r2,
           regionsIns=brachy_insertions_r2,
           random=FALSE)

simulateSV(output="./brachy/r3",
           genome=brachy,
           regionsDels=brachy_deletions_r3,
           regionsIns=brachy_insertions_r3,
           random=FALSE)

simulateSV(output="./brachy/r4",
           genome=brachy,
           regionsDels=brachy_deletions_r4,
           regionsIns=brachy_insertions_r4,
           random=FALSE)

simulateSV(output="./brachy/r5",
           genome=brachy,
           regionsDels=brachy_deletions_r5,
           regionsIns=brachy_insertions_r5,
           random=FALSE)
