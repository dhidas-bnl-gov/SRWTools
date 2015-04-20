#!/usr/bin/env python


import sys
import re
import matplotlib.pyplot as plt

for InFileName in sys.argv[1:]:
  print InFileName

  f = open(InFileName, 'r')
  
  m = re.search('X(.*)_Y(.*)_gap(.*).txt', InFileName)
  
  AtX   = m.group(1)
  AtY   = m.group(2)
  AtGap = m.group(3)
  
  
  Range_Z = [-850, 850]
  
  Z  = []
  Bx = []
  By = []
  Bz = []
  
  for l in f:
    v = l.split()
    Z.append(v[0])
    Bx.append(v[1])
    By.append(v[2])
    Bz.append(v[3])
  
  
  
  plt.figure()
  plt.subplot(311)
  plt.plot(Z, Bx)
  plt.xlim(Range_Z)
  plt.ylabel(r'$B_x (T)$')
  plt.xlabel("Z (mm)")
  
  plt.subplot(312)
  plt.plot(Z, By)
  plt.xlim(Range_Z)
  plt.ylabel(r'$B_y (T)$')
  plt.xlabel("Z (mm)")
  
  plt.subplot(313)
  plt.plot(Z, Bz)
  plt.xlim(Range_Z)
  plt.ylabel(r'$B_z (T)$')
  plt.xlabel("Z (mm)")
  
  
  OutFileName = 'FieldPlots_X' + AtX + '_Y' + AtY + '_gap' + AtGap + '.png'
  plt.savefig(OutFileName)
  
  
