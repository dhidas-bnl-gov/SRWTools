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


Bx2 = []
By2 = []
for i in range( len(Z) ):
  Bx2.append( IntegralBx[i] )
  By2.append( IntegralBy[i] )

Bx2Integral = IntegralVector(Z, Bx2)
By2Integral = IntegralVector(Z, By2)

print Bx2Integral[-1:][0] * 1e4, By2Integral[-1:][0] * 1e4, '  Gm^2'

plt.subplot(223)
plt.plot(Z, Bx2Integral)
plt.subplot(224)
plt.plot(Z, By2Integral)
plt.show()


#tck = interpolate.splrep(Z, By, s=0)
#xnew = numpy.arange(Z[0], Z[-1:][0], (Z[-1:][0] - Z[0]) / 100000.)
#ynew = interpolate.splev(xnew, tck, der=0)

#plt.plot(xnew, ynew)
#plt.show()



#Integralynew = IntegralVector(xnew, ynew)

#print Integralynew[-1:][0]
