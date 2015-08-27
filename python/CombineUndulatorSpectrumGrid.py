#!/usr/bin/env python

import sys


N = len(sys.argv) - 1
print 'N files:', N



Seeds = []
Values = []

for i in range(1, N + 1):
  FileName = sys.argv[i]
  print 'Adding:', FileName

  fi = open(FileName, 'r')

  Seed = int(fi.readline())

  # Check for seeing this seed
  if Seed in Seeds:
    print 'WARNING: Seed is already in the seed list.  skipping file:', FileName
    break

  j = 0
  for l in fi:
    l = l.strip()
    [x, y] = map(float, l.split())

    if i == 1:
      Values.append([x, y / N])
    else:
      Values[j][1] += y / N

    j += 1


#BaseFileName = os.path.splitext(os.path.basename(InFileName))[0]
#OutFileName = BaseFileName + '_Average.dat'
OutFileName = 'Average.dat'
fo = open(OutFileName, 'w')
fo.write('0\n')
for v in Values:
  fo.write(str(v[0]) + ' ' + str(v[1]) + '\n')
fo.close()
