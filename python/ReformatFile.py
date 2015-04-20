#!/usr/bin/env python

import sys

fi = open(sys.argv[1], 'r')
fo = open(sys.argv[2], 'w')


for l in fi:
  if float( l.split()[0] ) > -850 and float( l.split()[0] ) < 850:
    fo.write(' '.join(l.split()[1:]) + '\n')
