sed -i.bak 1d "insertions_${fname}.bed"
awk 'BEGIN {FS=OFS="\t"} {print "chr"$1,$2,$3,"chr"$7,$8,$9}' "insertions_${fname}.bed" > "${circos_data}_${fname}.txt"
