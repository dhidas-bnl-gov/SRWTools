#!/usr/bin/env python

import numpy
from array import array
import sys
from SRWToolsUtil import *

from ROOT import TLine

# Input and output files are given on command line
InFileName = sys.argv[1]
OutFileName = sys.argv[2]

# Open the input and output files.  Output file is ROOT format
fi = open(InFileName, 'r')
fo = TFile(OutFileName, 'recreate')

# Create a TTree for data 
t  = TTree('mmlabdata', 'mmlabdata')

# Create variables for TTree filling
tZ  = numpy.zeros(1, dtype=float)
tBx = numpy.zeros(1, dtype=float)
tBy = numpy.zeros(1, dtype=float)
tBz = numpy.zeros(1, dtype=float)

# Create branches for the tree
t.Branch('Z',  tZ,  'Z')
t.Branch('Bx', tBx, 'Bx')
t.Branch('By', tBy, 'By')
t.Branch('Bz', tBz, 'Bz')

# Create arrays for program use
Z  = array('d')
Bx = array('d')
By = array('d')
Bz = array('d')


# Loop over all entries in the input file
for l in fi:
  [tZ[0], tBx[0], tBy[0], tBz[0]] = map(float, l.split())

  # Make in SI unites and fill the Tree
  tZ[0] /= 1000.
  t.Fill()

  # Fill program use arrays with data
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

# TGraphs for measured Bx By and Bz vs Z
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

# Write graphs to open TFile
gBx.Write()
gBy.Write()
gBz.Write()


# Copy list and pop off each side of MaxBy, Average max's.
# Here I am avoiding the fringe on both sides
MaxByZChopped = MaxByZ[:]
MaxByChopped  = MaxBy[:]
for i in range (6):
  MaxByChopped.pop()
  MaxByZChopped.pop()
  MaxByChopped.pop(0)
  MaxByZChopped.pop(0)

# just save absolute value
MaxByChopped = map(abs, MaxByChopped)

# For comparing the number of peaks found and taken off on the sides
print 'MaxBy MaxByChopped:', len(MaxByChopped), len(MaxBy)

# Number of periods based on measured data
NPERIODS = (len(MaxBy) - 4) / 2 

# Calculate the agerave period from Chopped
HalfPeriodList = []
for i in range(1, len(MaxByZ)):
  HalfPeriodList.append( MaxByZ[i] - MaxByZ[i-1])

# The period length as seen in data
PERIOD_LENGTH = numpy.mean(HalfPeriodList) * 2
print 'Period length measured: ', PERIOD_LENGTH, '[m]'


# Average max field
AvgMaxBy = numpy.mean(MaxByChopped)
StdMaxBy = numpy.std(MaxByChopped)
print 'Max By Average: ', AvgMaxBy, ' +/- ', StdMaxBy, ' [T]'

# Magnetic center point of undulator
UNDULATOR_ZCENTER = numpy.sum(MaxByZChopped)
print 'UNDULATOR_ZCENTER: ', UNDULATOR_ZCENTER, '[m]'

# Start and end of Undulator in Z
UNDULATOR_ZSTART = UNDULATOR_ZCENTER - ( PERIOD_LENGTH * (NPERIODS + 4) / 2)
UNDULATOR_ZEND   = UNDULATOR_ZCENTER + ( PERIOD_LENGTH * (NPERIODS + 4) / 2)
UNDULATOR_LENGTH = UNDULATOR_ZEND - UNDULATOR_ZSTART
print 'UNDULATOR ZSTART ZEND LENGTH:', UNDULATOR_ZSTART, UNDULATOR_ZEND, UNDULATOR_LENGTH

lStart = TLine(UNDULATOR_ZSTART, -1, UNDULATOR_ZSTART, 1)
lStart.SetLineColor(2)
lCenter = TLine(UNDULATOR_ZCENTER, -1, UNDULATOR_ZCENTER, 1)
lCenter.SetLineColor(2)
lEnd = TLine(UNDULATOR_ZEND, -1, UNDULATOR_ZEND, 1)
lEnd.SetLineColor(2)

cBxLines = TCanvas('Undulator_Bx')
gBx.Draw("ACP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cBxLines.Write()

cByLines = TCanvas('Undulator_By')
gBy.Draw("ACP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cByLines.Write()

cBzLines = TCanvas('Undulator_Bz')
gBz.Draw("ACP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cBzLines.Write()


# Peak By to use in calculations
PEAK_BY = AvgMaxBy

# Graph the distribution and fit a line while we're at it
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
magFldCnt_Data.arMagFld[0] = ReadHallProbeDataSRW(InFileName, UNDULATOR_ZSTART, UNDULATOR_ZEND)

# Field interpolation method
magFldCnt_Data.arMagFld[0].interp = 4

# ID center
magFldCnt_Data.arXc[0] = 0.0
magFldCnt_Data.arYc[0] = UNDULATOR_ZCENTER


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
zcID = UNDULATOR_ZCENTER


und = SRWLMagFldU([SRWLMagFldH(1, 'v', undBy, phBy, sBy, 1), SRWLMagFldH(1, 'h', undBx, phBx, sBx, 1)], undPer, numPer)
magFldCnt_Ideal = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements

# Get the electron trajectory
partTraj_Ideal = GetElectronTrajectory(magFldCnt_Ideal, -1, 1)

# Get the Z Values for the plot
ZValues = [float(x) * ((partTraj_Ideal.ctEnd - partTraj_Ideal.ctStart) / float(partTraj_Ideal.np)) for x in range(0, partTraj_Ideal.np)]


gElectronX_Ideal = TGraph( len(ZValues), array('d', ZValues), partTraj_Ideal.arX)
gElectronY_Ideal = TGraph( len(ZValues), array('d', ZValues), partTraj_Ideal.arY)

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


exit(0)


def dimdelta (X) :
  "step size based on last minus first divide by N-1"

  N = len(X) - 1
  return (X[-1] - X[0]) / N


def dimoffset (X) :
  "return the 0th element.  Should be Z im m"
  return X[0]


# Field integrals
FInt1Bx = IntegralVector(Z, Bx)
FInt2Bx = IntegralVector(Z, FInt1Bx)
FInt1By = IntegralVector(Z, By)
FInt2By = IntegralVector(Z, FInt1By)

auxI1X = FInt1Bx[-1]
auxI2X = FInt2Bx[-1]
auxI1Y = FInt1By[-1]
auxI2Y = FInt2By[-1]



# Get sRangeY = length of undulator, not entire field
npAuxFldInt = len(Bx)
sRangeX = (npAuxFldInt - 1) * dimdelta(Z)

npAuxFldInt = len(By)
sRangeY = (npAuxFldInt - 1) * dimdelta(Z)




# Distance between kicks
DistanceBetweenKicks = PERIOD_LENGTH * NPERIODS
sCenX = 0
sCenY = 0

# RMS Kick Length
rmsLenKick = 0.0001

KickEntryHorizontal =  0.5 * (sRangeY / DistanceBetweenKicks - 1) * auxI1Y - auxI2Y / DistanceBetweenKicks
KickExitHorizontal  = -0.5 * (sRangeY / DistanceBetweenKicks + 1) * auxI1Y + auxI2Y / DistanceBetweenKicks
KickEntryVertical   =  0.5 * (sRangeX / DistanceBetweenKicks - 1) * auxI1X - auxI2X / DistanceBetweenKicks
KickExitVertical    = -0.5 * (sRangeX / DistanceBetweenKicks + 1) * auxI1X + auxI2X / DistanceBetweenKicks


print 'DistanceBetweenKicks: ', DistanceBetweenKicks
print 'Kick EnH ExH EnV ExV: ', KickEntryHorizontal, KickExitHorizontal, KickEntryVertical, KickExitVertical

def AddToField (X, B, s0, sigma, intgr) :
  "Add gaussian to field"

  for i in range(len(X)):
    B[i] += sRGssn( X[i], s0, sigma, intgr )

def sRGssn(s, s0, sigma, intgr):
  "gaussian"

  t=(s-s0)/sigma
  return (intgr*0.3989422804/sigma)*exp(-0.5*t*t)


HalfDistBwKicks = 0.5 * DistanceBetweenKicks

AddToField(Z, By, sCenY - HalfDistBwKicks, rmsLenKick, KickEntryHorizontal)
AddToField(Z, By, sCenY + HalfDistBwKicks, rmsLenKick, KickExitHorizontal)
AddToField(Z, Bx, sCenX - HalfDistBwKicks, rmsLenKick, KickEntryVertical)
AddToField(Z, Bx, sCenX + HalfDistBwKicks, rmsLenKick, KickExitVertical)


#AddToField(Z, By, -sRangeY/2 + 3 * rmsLenKick, rmsLenKick, KickEntryHorizontal)
#AddToField(Z, By,  sRangeY/2 - 3 * rmsLenKick, rmsLenKick, KickExitHorizontal)
#AddToField(Z, Bx, -sRangeX/2 + 3 * rmsLenKick, rmsLenKick, KickEntryVertical)
#AddToField(Z, Bx,  sRangeX/2 - 3 * rmsLenKick, rmsLenKick, KickExitVertical)


exit(0)





magFldCnt_Corr = SRWLMagFldC()
magFldCnt_Corr.allocate(1)


# Make magnetic field object for corrected field
magFldCnt_Corr.arMagFld[0] = SRWLMagFld3D( array('d', Bx), array('d', By), array('d', Bz), 1, 1, len(Z), 0.0, 0.0, Z[-1] - Z[0], 1, 1, None, None, _arZ=array('d', Z))

# Field interpolation method
magFldCnt_Corr.arMagFld[0].interp = 4

# ID center
magFldCnt_Corr.arXc[0] = 0.0
magFldCnt_Corr.arYc[0] = 0.0


# Number of reps of field
magFldCnt_Corr.arMagFld[0].nRep = 1

# Center in Z of ID
magFldCnt_Corr.arZc[0] = 0.0








# Get the electron trajectory
partTraj_Corr = GetElectronTrajectory(magFldCnt_Corr, -1, 1)







# Get the Z Values for the plot
ZValues_Corr = [float(x) * ((partTraj_Corr.ctEnd - partTraj_Corr.ctStart) / float(partTraj_Corr.np)) for x in range(0, partTraj_Corr.np)]


gElectronX_Corr = TGraph( len(ZValues_Corr), array('d', ZValues_Corr), partTraj_Corr.arX)
gElectronY_Corr = TGraph( len(ZValues_Corr), array('d', ZValues_Corr), partTraj_Corr.arY)

gElectronX_Corr.SetTitle('Electron Trajectory in X')
gElectronX_Corr.GetXaxis().SetTitle('Z Position [m]')
gElectronX_Corr.GetYaxis().SetTitle('X Position [m]')
gElectronX_Corr.SetName("ElectronX_Corr")

gElectronY_Corr.SetTitle('Electron Trajectory in Y')
gElectronY_Corr.GetXaxis().SetTitle('Z Position [m]')
gElectronY_Corr.GetYaxis().SetTitle('Y Position [m]')
gElectronY_Corr.SetName("ElectronY_Corr")

gElectronX_Corr.Write()
gElectronY_Corr.Write()




[SpectrumDataX, SpectrumDataY] = GetUndulatorSpectrum(magFldCnt_Data)
gSpectrumData = TGraph( len(SpectrumDataX), array('d', SpectrumDataX), array('d', SpectrumDataY) )
gSpectrumData.SetName('SpectrumData')
gSpectrumData.SetTitle('Simulated B-field Spectrum from Field Measurements')
gSpectrumData.GetXaxis().SetTitle('Energy [eV]')
gSpectrumData.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumData.Write()


[SpectrumCorrX, SpectrumCorrY] = GetUndulatorSpectrum(magFldCnt_Corr)
gSpectrumCorr = TGraph( len(SpectrumCorrX), array('d', SpectrumCorrX), array('d', SpectrumCorrY) )
gSpectrumCorr.SetName('SpectrumCorr')
gSpectrumCorr.SetTitle('Simulated B-field Spectrum from Field Measurements')
gSpectrumCorr.GetXaxis().SetTitle('Energy [eV]')
gSpectrumCorr.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumCorr.Write()



fo.Write()
fo.Close()

