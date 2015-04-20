#!/usr/bin/env python


import sys
import re
import matplotlib.pyplot as plt

InFileName = sys.argv[1]

f = open(InFileName, 'r')


Range_Z = [-850, 850]

Z  = []
Bx = []
By = []
Bz = []

for l in f:
  v = l.split()
  if ( float(v[0]) > 800.0 and float(v[0]) < 850.0):
    Z.append(float(v[0]) * 1e3)


dZ = []
X = []
for i in range(len(Z)):
  if i > 1:
    dZ.append(Z[i] - Z[i-1])
    X.append(i-1)


plt.figure()
plt.subplot(211)
plt.plot(X, dZ)
plt.ylabel('dZ (um)')
plt.xlabel("move")

plt.subplot(212)
plt.hist(dZ)



OutFileName = 'dZ.png'
plt.savefig(OutFileName)
  
  
