from __future__ import division
import sys
import math

lengths = []
vari = 0
for line in sys.stdin:
    if line.startswith('@'):
        pass
    else:
        line = line.rsplit()
        tlen = int(line[8])
        if tlen > 0:
            lengths.append(tlen)
            vari += (tlen*tlen)
        else:
            pass

N = len(lengths)
S = sum(lengths)
mean = int(S/N)
val = ((vari-(mean*mean*N))/(N-1))
std = int(math.sqrt(val))

print mean, std
