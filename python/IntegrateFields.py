#!/usr/bin/env python

import sys
from SRWToolsUtil import *
import matplotlib.pyplot as plt

import numpy
from scipy import interpolate


# grab data from file
[Z, Bx, By, Bz] = ReadHallProbeData(sys.argv[1], -1, 1)


# trap rule for first integral
IntegralBx = IntegralVector(Z, Bx)
IntegralBy = IntegralVector(Z, By)

plt.subplot(221)
plt.plot(Z, IntegralBx)
plt.subplot(222)
plt.plot(Z, IntegralBy)

print IntegralBx[-1:][0] * 1e4, IntegralBy[-1:][0] * 1e4, '  Gm'



Integral2Bx = IntegralVector(Z, IntegralBx)
Integral2By = IntegralVector(Z, IntegralBy)

print Bx2Integral[-1:][0] * 1e4, By2Integral[-1:][0] * 1e4, '  Gm^2'

plt.subplot(223)
plt.plot(Z, Bx2Integral)
plt.subplot(224)
plt.plot(Z, By2Integral)
plt.show()



