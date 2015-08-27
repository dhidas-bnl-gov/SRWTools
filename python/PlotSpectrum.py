#!/usr/bin/env python

import sys
from ROOT import TGraph, TCanvas


fi = open(sys.argv[1])

fi.readline()

g = TGraph()
g.SetTitle('Undulator Spectrum')

n = 0
for l in fi:
  n += 1
  [x, y] = map(float, l.split())
  g.Set(n)
  g.SetPoint(n-1, x, y)


c = TCanvas()
c.cd()
g.Draw('ACP')
c.SaveAs('Average.pdf')
