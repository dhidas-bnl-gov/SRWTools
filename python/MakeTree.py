#!/usr/bin/env python

import numpy
from array import array
import sys
from SRWToolsUtil import *

InFileName = sys.argv[1]
OutFileName = sys.argv[2]

fi = open(InFileName, 'r')
fo = TFile(OutFileName, 'recreate')
t  = TTree('mmlabdata', 'mmlabdata')

tZ  = numpy.zeros(1, dtype=float)
tBx = numpy.zeros(1, dtype=float)
tBy = numpy.zeros(1, dtype=float)
tBz = numpy.zeros(1, dtype=float)

t.Branch('Z',  tZ,  'Z')
t.Branch('Bx', tBx, 'Bx')
t.Branch('By', tBy, 'By')
t.Branch('Bz', tBz, 'Bz')

Z  = array('d')
Bx = array('d')
By = array('d')
Bz = array('d')


for l in fi:
  [tZ[0], tBx[0], tBy[0], tBz[0]] = map(float, l.split())

  # make in SI unites
  tZ[0] /= 1000.
  t.Fill()

  Z.append(tZ[0])
  Bx.append(tBx[0])
  By.append(tBy[0])
  Bz.append(tBz[0])


# Get the max and mins in the field
[MaxListInd, MaxListBy] = FindMaxAndMins(Z, By)
[MaxByZ, MaxBy] = FindMinMaxFromFit(MaxListInd, Z, By)

# Save a graph of the points found
gMaxBy = TGraph(len(MaxByZ), array('d', MaxByZ), array('d', MaxBy))
gMaxBy.SetName('MaxBy')
gMaxBy.SetTitle('Calculated Max By')
gMaxBy.GetXaxis().SetTitle('Position [m]')
gMaxBy.GetYaxis().SetTitle('B_{Y} [T]')
gMaxBy.Write()

gBx = TGraph( len(Z), Z, Bx )
gBy = TGraph( len(Z), Z, By )
gBz = TGraph( len(Z), Z, Bz )

gBx.SetTitle('Measured magnetic field B_{X}')
gBx.GetXaxis().SetTitle('Position [m]')
gBx.GetYaxis().SetTitle('#B_{x} [T]')
gBx.SetName("Bx")

gBy.SetTitle('Measured magnetic field B_{Y}')
gBy.GetXaxis().SetTitle('Position [m]')
gBy.GetYaxis().SetTitle('B_{Y} [T]')
gBy.SetMarkerStyle(26)
gBy.SetName("By")

gBz.SetTitle('Measured magnetic field B_{Z}')
gBz.GetXaxis().SetTitle('Position [m]')
gBz.GetYaxis().SetTitle('B_{Z} [T]')
gBz.SetName("Bz")

gBx.Write()
gBy.Write()
gBz.Write()


# copy list and pop 4 off each side of MaxBy, Average max's
MaxByZChopped = MaxByZ[:]
MaxByChopped  = MaxBy[:]

for i in range (6):
  MaxByChopped.pop()
  MaxByZChopped.pop()
  MaxByChopped.pop(0)
  MaxByZChopped.pop(0)

# just save absolute value
MaxByChopped = map(abs, MaxByChopped)

print len(MaxByChopped), len(MaxByZChopped), len(MaxBy), len(MaxByZ)

# number of periods based on measured data
NPERIODS = (len(MaxBy) - 4) / 2

# Calculate the agerave period from Chopped
HalfPeriodList = []
for i in range(1, len(MaxByZ)):
  HalfPeriodList.append( MaxByZ[i] - MaxByZ[i-1])

# The period length as seen in data
PERIOD_LENGTH = numpy.mean(HalfPeriodList) * 2
print 'Period length measured: ', PERIOD_LENGTH


# average max field
AvgMaxBy = numpy.mean(MaxByChopped)
StdMaxBy = numpy.std(MaxByChopped)
print 'Max By Average: ', AvgMaxBy, ' +/- ', StdMaxBy, ' [T]'

# peak By to use in calculations
PEAK_BY = AvgMaxBy

# graph the distribution and fit a line while we're at it
gMaxByChopped = TGraph(len(MaxByZChopped), array('d', MaxByZChopped), array('d', MaxByChopped))
gMaxByChopped.SetName('MaxByChopped')
gMaxByChopped.SetTitle('Calculated Max By')
gMaxByChopped.GetXaxis().SetTitle('Position [m]')
gMaxByChopped.GetYaxis().SetTitle('B_{Y} [T]')
fline = TF1('fline', 'pol1', -3000, 3000)
gMaxByChopped.Fit(fline, 'q')
gMaxByChopped.SetMarkerStyle(32)
gMaxByChopped.Write()
print 'Max By linear fit parameters [0]+[1]x::', fline.GetParameter(0), fline.GetParameter(1)










# Defining Magnetic Field container:
magFldCnt_Data = SRWLMagFldC()
magFldCnt_Data.allocate(1)


# Read data from file and make mag field object
magFldCnt_Data.arMagFld[0] = ReadHallProbeDataSRW(InFileName)

# Field interpolation method
magFldCnt_Data.arMagFld[0].interp = 4

# ID center
magFldCnt_Data.arXc[0] = 0.0
magFldCnt_Data.arYc[0] = 0.0


# Number of reps of field
magFldCnt_Data.arMagFld[0].nRep = 1

# Center in Z of ID
magFldCnt_Data.arZc[0] = 0.0

# Get the electron trajectory
partTraj = GetElectronTrajectory(magFldCnt_Data, -1, 1)

# Get the Z Values for the plot
ZValues = [float(x) * ((partTraj.ctEnd - partTraj.ctStart) / float(partTraj.np)) for x in range(0, partTraj.np)]


gElectronX_Data = TGraph( len(ZValues), array('d', ZValues), partTraj.arX)
gElectronY_Data = TGraph( len(ZValues), array('d', ZValues), partTraj.arY)

gElectronX_Data.SetTitle('Electron Trajectory in X')
gElectronX_Data.GetXaxis().SetTitle('Z Position [m]')
gElectronX_Data.GetYaxis().SetTitle('X Position [m]')
gElectronX_Data.SetName("ElectronX_Data")

gElectronY_Data.SetTitle('Electron Trajectory in Y')
gElectronY_Data.GetXaxis().SetTitle('Z Position [m]')
gElectronY_Data.GetYaxis().SetTitle('Y Position [m]')
gElectronY_Data.SetName("ElectronY_Data")

gElectronX_Data.Write()
gElectronY_Data.Write()








# Defining Magnetic Field container:

# define the undulator specs
numPer = NPERIODS
undPer = PERIOD_LENGTH
undBx = 0.0
undBy = PEAK_BY
phBx = 0
phBy = 0
sBx = 1
sBy = 1
xcID = 0
ycID = 0
zcID = 0


und = SRWLMagFldU([SRWLMagFldH(1, 'v', undBy, phBy, sBy, 1), SRWLMagFldH(1, 'h', undBx, phBx, sBx, 1)], undPer, numPer)
magFldCnt_Ideal = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements

# Get the electron trajectory
partTraj = GetElectronTrajectory(magFldCnt_Ideal, -1, 1)

# Get the Z Values for the plot
ZValues = [float(x) * ((partTraj.ctEnd - partTraj.ctStart) / float(partTraj.np)) for x in range(0, partTraj.np)]


gElectronX_Ideal = TGraph( len(ZValues), array('d', ZValues), partTraj.arX)
gElectronY_Ideal = TGraph( len(ZValues), array('d', ZValues), partTraj.arY)

gElectronX_Ideal.SetTitle('Electron Trajectory in X')
gElectronX_Ideal.GetXaxis().SetTitle('Z Position [m]')
gElectronX_Ideal.GetYaxis().SetTitle('X Position [m]')
gElectronX_Ideal.SetName("ElectronX_Ideal")

gElectronY_Ideal.SetTitle('Electron Trajectory in Y')
gElectronY_Ideal.GetXaxis().SetTitle('Z Position [m]')
gElectronY_Ideal.GetYaxis().SetTitle('Y Position [m]')
gElectronY_Ideal.SetName("ElectronY_Ideal")

gElectronX_Ideal.Write()
gElectronY_Ideal.Write()





[SpectrumIdealX, SpectrumIdealY] = GetUndulatorSpectrum(magFldCnt_Ideal)
gSpectrumIdeal = TGraph( len(SpectrumIdealX), array('d', SpectrumIdealX), array('d', SpectrumIdealY) )
gSpectrumIdeal.SetName('SpectrumIdeal')
gSpectrumIdeal.SetTitle('Simulated B-field Spectrum')
gSpectrumIdeal.GetXaxis().SetTitle('Energy [eV]')
gSpectrumIdeal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumIdeal.Write()





# First field integral
FInt1Bx = IntegralVector(Z, Bx)
print FInt1Bx[-1]
exit(0)



[SpectrumX, SpectrumY] = GetUndulatorSpectrum(magFldCnt_Data)
gSpectrum = TGraph( len(SpectrumX), array('d', SpectrumX), array('d', SpectrumY) )
gSpectrum.SetName('SpectrumData')
gSpectrum.SetTitle('Simulated B-field Spectrum from Field Measurements')
gSpectrum.GetXaxis().SetTitle('Energy [eV]')
gSpectrum.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrum.Write()


fo.Write()
fo.Close()

