#!/usr/bin/env python

import sys
from SRWToolsUtil import *
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d


# grab data from file
[Z, Bx, By, Bz] = ReadHallProbeData(sys.argv[1])


# trap rule for first integral
IntegralBx = [0.0]
IntegralBy = [0.0]
for i in range(1, len(Z)):
  dz = Z[i] - Z[i-1]
  trap = dz * (Bx[i-1] + 0.5 * (By[i] - By[i-1]))
  IntegralBx.append(IntegralBx[-1:][0] + trap)
  trap = dz * (By[i-1] + 0.5 * (By[i] - By[i-1]))
  IntegralBy.append(IntegralBy[-1:][0] + trap)


#plt.plot(Z, IntegralBx)
#plt.show()


#plt.plot(Z, IntegralBy)
#plt.show()

print IntegralBx[-1:][0], IntegralBy[-1:][0]





