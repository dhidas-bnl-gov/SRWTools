#!/usr/bin/env python

import sys
from SRWToolsUtil import *
import matplotlib.pyplot as plt

import numpy
from scipy import interpolate


# grab data from file
[Z, Bx, By, Bz] = ReadHallProbeData(sys.argv[1])

fo = open(sys.argv[2], 'w')



ZNew = numpy.arange(Z[0], Z[-1:][0], (Z[-1:][0] - Z[0]) / 100000.)

tckX = interpolate.splrep(Z, Bx, s=0)
BxNew = interpolate.splev(ZNew, tckX, der=0)

tckY = interpolate.splrep(Z, By, s=0)
ByNew = interpolate.splev(ZNew, tckY, der=0)

tckZ = interpolate.splrep(Z, Bz, s=0)
BzNew = interpolate.splev(ZNew, tckZ, der=0)


for i in range( len(ZNew) ):
  l = map(str, [ZNew[i], BxNew[i], ByNew[i], BzNew[i]])
  fo.write(' '.join(l) + '\n')
