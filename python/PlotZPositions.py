#!/usr/bin/env python

from SRWToolsUtil import *

import sys
import os.path
import re
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
  print 'Usage: ', sys.argv[0], ' [InFile] [OutFile]'
  sys.exit()

InFileName = sys.argv[1]
OutFileName = sys.argv[2]

fi = open(InFileName, 'r')

if os.path.exists(OutFileName):
  print 'File exists already: ', OutFileName
  sys.exit()

fo = open(OutFileName, 'w')



Range_Z = [-850., 850.]

Z  = []
Bx = []
By = []
Bz = []

for l in fi:
  v = l.split()
  if ( float(v[0]) > Range_Z[0] and float(v[0]) < Range_Z[1]):
    Z.append(float(v[0]) * 1e3)


dZ = []
X = []
for i in range(len(Z)):
  if i > 1:
    dZ.append(Z[i] - Z[i-1])
    X.append(Z[i] / 1000.)


plt.figure()
plt.subplot(211)
plt.plot(X, dZ)
plt.ylabel('dZ (um)')
plt.xlabel("move")

plt.subplot(212)
plt.hist(dZ)



plt.savefig(OutFileName)
  
  
