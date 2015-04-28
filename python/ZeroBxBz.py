#!/usr/bin/env python

import sys
from SRWToolsUtil import *



# grab data from file
[Z, Bx, By, Bz] = ReadHallProbeData(sys.argv[1])

fo = open(sys.argv[2], 'w')


for i in range( len(Z) ):
  l = map(str, [Z[i], 0, By[i], 0])
  fo.write(' '.join(l) + '\n')
