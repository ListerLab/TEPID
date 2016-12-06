# Documentation for scripts

```
merge_insertions.py
```

Merges each TE insertion in a population if the coordinates are within +/- 100 bp and the TE name is the same and creates a list of accessions containing that insertion.

input:  
  `-f` Filename prefix for the files to be merged. Will likely be one of `insertions`, `ambiguous_insertions`, or `refined_insertions`
  
output:  
  * File with two sets of coordinates for each variant (origin and insertion coordinates), along with a list of all accessions containing the insertion.  
  * File with one set of coordinates (the insertion coordinates), of the same format as the TE deletion identification output, to allow the two datasets to be merged.


```
merge_deletions.py
```

Merges TE deletions within a population and creates a list of accessions containing the TE deletion.

input:  
  `-f` Filename prefix for the files to be merged. Will likely be one of `deletions`, `ambiguous_deletions` or `refined_deletions`
  
output:  
  * A single file withn the coordinates of each TE deletion and a list of accessions containing that deletion.
  
```
flip_deletions.py
```

Inverts the list of accessions containing a TE deletion, so that the list is the accessions that contain a TE present at the given coordinates, in order to be consistent with the data from the TE insertion discovery steps.

input:  
  `-s` list of sample names (text file with one accession per line)  
  `-d` file containing the merged deletions  
  `-r` name of the reference accession  
  `-o` output filename
  
output:
  * A file formatted the same way as the output of `merge_insertions.py`, allowing the two to be concatenated.
  
```
genotype.py
```

Only to be run after the refinement step. Creates a list of accessions that contain the TE variant, and another list that does not contain the TE variant. There may be accessions that are not in either list, where there was insufficient coverage at that site to make a call.

input:  
  `-d` run on deletions  
  `-i` run on insertions  
  `-a` ambiguous variants filename  
  `-r` name of the reference accession  
  
output:  
  * prints to stdout  
  * similar file format to the output of `merge_insertions.py` but has a additional list of accessions that do not contain each variant.
