#!/usr/bin/env python

import os
import numpy
from array import array
import sys

sys.argv.append('-b-')

from SRWToolsUtil import *

# extra root imports
from ROOT import TLine, TH1F, TH2F, TGraph, TCanvas, TFile








# Input and output files are given on command line
InFileName = sys.argv[1]
print 'Input file:', InFileName


# Name of output files for this section
BaseFileName = os.path.splitext(os.path.basename(InFileName))[0]
OutFileNameRoot = BaseFileName + '_MMFAnalysis.root'
print 'Output file for MMFAnalysis: ', OutFileNameRoot


# Open the input and output ROOT file.
fIN = open(InFileName, 'r')
fROOT = TFile(OutFileNameRoot, 'recreate')
fROOT.cd()


# Create arrays for magnetic field.  These contain the field data for the entire input file
Z  = array('d')
Bx = array('d')
By = array('d')
Bz = array('d')


# Loop over all entries in the input file
for l in fIN:
  Data = map(float, l.split())

  # Fill program use arrays with data
  Z.append(Data[0] / 1000.)
  Bx.append(Data[1])
  By.append(Data[2])
  Bz.append(Data[3])

# Close file
fIN.close()


# Get the max and mins in the field
[MaxListInd, MaxListBy] = FindMaxAndMins(Z, By)
[MaxByZ, MaxBy] = FindMinMaxFromFit(MaxListInd, Z, By, None)


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

# Just save absolute value
MaxByChoppedAbs = map(abs, MaxByChopped)

# For comparing the number of peaks found and taken off on the sides
print 'MaxByChopped MaxBy:', len(MaxByChopped), len(MaxBy)

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
AvgMaxBy = numpy.mean(MaxByChoppedAbs)
StdMaxBy = numpy.std(MaxByChoppedAbs)
print 'Max By Average: ', AvgMaxBy, ' +/- ', StdMaxBy, ' [T]'

# Magnetic center point of undulator
UNDULATOR_ZCENTER = numpy.sum(MaxByZChopped)
print 'UNDULATOR_ZCENTER: ', UNDULATOR_ZCENTER, '[m]'

# Start and end of Undulator in Z
UNDULATOR_ZSTART = UNDULATOR_ZCENTER - ( PERIOD_LENGTH * (NPERIODS + 4) / 2)
UNDULATOR_ZEND   = UNDULATOR_ZCENTER + ( PERIOD_LENGTH * (NPERIODS + 4) / 2)
UNDULATOR_LENGTH = UNDULATOR_ZEND - UNDULATOR_ZSTART
print 'UNDULATOR ZSTART ZEND LENGTH:', UNDULATOR_ZSTART, UNDULATOR_ZEND, UNDULATOR_LENGTH

# Peak By to use in calculations
PEAK_BY = AvgMaxBy



# Undulator fields and START, CENTER, END lines
lStart = TLine(UNDULATOR_ZSTART, -1, UNDULATOR_ZSTART, 1)
lStart.SetLineColor(2)
lCenter = TLine(UNDULATOR_ZCENTER, -1, UNDULATOR_ZCENTER, 1)
lCenter.SetLineColor(2)
lEnd = TLine(UNDULATOR_ZEND, -1, UNDULATOR_ZEND, 1)
lEnd.SetLineColor(2)

# Canvas for Bx
cBxLines = TCanvas('Undulator_Bx')
gBx.Draw("ALP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cBxLines.Write()

# Canvas for By
cByLines = TCanvas('Undulator_By')
gBy.SetMarkerStyle(3)
gBy.SetMarkerSize(0.3)
gBy.Draw("ALP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cByLines.Write()

# Canvas for Bz
cBzLines = TCanvas('Undulator_Bz')
gBz.Draw("ALP")
lStart.Draw()
lCenter.Draw()
lEnd.Draw()
cBzLines.Write()


# Graph the distribution and fit a line while we're at it
gMaxByChoppedAbs = TGraph(len(MaxByZChopped), array('d', MaxByZChopped), array('d', MaxByChoppedAbs))
gMaxByChoppedAbs.SetName('MaxByChoppedAbs')
gMaxByChoppedAbs.SetTitle('Calculated Max By')
gMaxByChoppedAbs.GetXaxis().SetTitle('Position [m]')
gMaxByChoppedAbs.GetYaxis().SetTitle('B_{Y} [T]')
fline = TF1('fline', 'pol1', -3000, 3000)
gMaxByChoppedAbs.Fit(fline, 'q')
gMaxByChoppedAbs.SetMarkerStyle(32)
print 'Max By linear fit parameters [0]+[1]x::', fline.GetParameter(0), fline.GetParameter(1)

# Write to output root file
gMaxByChoppedAbs.Write()




# Trajectory start and stop values
TRAJECTORY_START = UNDULATOR_ZSTART - (UNDULATOR_LENGTH * 0.10)
TRAJECTORY_STOP  = UNDULATOR_ZEND   + (UNDULATOR_LENGTH * 0.10)




# Get undulator and Ideal magnetic field based off of average peak By field
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
PartTraj_Ideal = GetElectronTrajectory(magFldCnt_Ideal, TRAJECTORY_START, TRAJECTORY_STOP)

# Get the Z Values for the plot
ZValues = numpy.linspace(PartTraj_Ideal.ctStart, PartTraj_Ideal.ctEnd, PartTraj_Ideal.np)


# Graphs for trajectory from Ideal field
gElectronX_Ideal = TGraph( len(ZValues), array('d', ZValues), PartTraj_Ideal.arX)
gElectronY_Ideal = TGraph( len(ZValues), array('d', ZValues), PartTraj_Ideal.arY)

gElectronX_Ideal.SetTitle('Ideal Electron Trajectory in X')
gElectronX_Ideal.GetXaxis().SetTitle('Z Position [m]')
gElectronX_Ideal.GetYaxis().SetTitle('X Position [m]')
gElectronX_Ideal.SetName("ElectronX_Ideal")

gElectronY_Ideal.SetTitle('Ideal Electron Trajectory in Y')
gElectronY_Ideal.GetXaxis().SetTitle('Z Position [m]')
gElectronY_Ideal.GetYaxis().SetTitle('Y Position [m]')
gElectronY_Ideal.SetName("ElectronY_Ideal")

gElectronX_Ideal.Write()
gElectronY_Ideal.Write()










# Defining Magnetic Field container:
magFldCnt_Data = SRWLMagFldC()
magFldCnt_Data.allocate(1)


# Read data from file and make mag field object
#magFldCnt_Data.arMagFld[0] = ReadHallProbeDataSRW(InFileName, UNDULATOR_ZSTART, UNDULATOR_ZEND)
magFldCnt_Data.arMagFld[0] = SRWLMagFld3D( array('d', Bx), array('d', By), array('d', Bz), 1, 1, len(Z), 0.0, 0.0, Z[-1] - Z[0], 1, 1, None, None, _arZ=array('d', Z))

# Field interpolation method
magFldCnt_Data.arMagFld[0].interp = 4

# ID center
magFldCnt_Data.arXc[0] = 0.0
magFldCnt_Data.arYc[0] = 0.0
magFldCnt_Data.arZc[0] = UNDULATOR_ZCENTER

# Number of reps of field
magFldCnt_Data.arMagFld[0].nRep = 1

# Get the electron trajectory
PartTraj_Data = GetElectronTrajectory(magFldCnt_Data, TRAJECTORY_START, TRAJECTORY_STOP)


# Get the Z Values for the plot
ZValues = numpy.linspace(PartTraj_Data.ctStart, PartTraj_Data.ctEnd, PartTraj_Data.np)

# Graphs for trajectory from uncorrected data
gElectronX_Data = TGraph( len(ZValues), array('d', ZValues), PartTraj_Data.arX)
gElectronY_Data = TGraph( len(ZValues), array('d', ZValues), PartTraj_Data.arY)

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









# Field integral vectors
FInt1Bx = IntegralVector(Z, Bx)
FInt2Bx = IntegralVector(Z, FInt1Bx)
FInt1By = IntegralVector(Z, By)
FInt2By = IntegralVector(Z, FInt1By)

# Field Integrals
auxI1X = FInt1Bx[-1]
auxI2X = FInt2Bx[-1]
auxI1Y = FInt1By[-1]
auxI2Y = FInt2By[-1]


# Get sRange = length of undulator, not entire field
npAuxFldInt = len(Bx)
sRangeX = (npAuxFldInt - 1) * StepSize(Z)
npAuxFldInt = len(By)
sRangeY = (npAuxFldInt - 1) * StepSize(Z)


# Distance between kicks and center X Y
DistanceBetweenKicks = UNDULATOR_ZEND - UNDULATOR_ZSTART
sCenX = 0
sCenY = 0

# RMS Kick Length in [m]
rmsLenKick = 0.050

# Correction kicks
KickEntryHorizontal =  0.5 * (sRangeY / DistanceBetweenKicks - 1) * auxI1Y - auxI2Y / DistanceBetweenKicks
KickExitHorizontal  = -0.5 * (sRangeY / DistanceBetweenKicks + 1) * auxI1Y + auxI2Y / DistanceBetweenKicks
KickEntryVertical   =  0.5 * (sRangeX / DistanceBetweenKicks - 1) * auxI1X - auxI2X / DistanceBetweenKicks
KickExitVertical    = -0.5 * (sRangeX / DistanceBetweenKicks + 1) * auxI1X + auxI2X / DistanceBetweenKicks

print 'DistanceBetweenKicks: %+.6E' % DistanceBetweenKicks
print 'Kick EntryHor %+.6E ExitHor %+.6E EntryVert %+.6E ExitVert %+.6E' % (KickEntryHorizontal, KickExitHorizontal, KickEntryVertical, KickExitVertical)


# Add kicks to magnetic field
AddToField(Z, By, UNDULATOR_ZCENTER - 0.5 * DistanceBetweenKicks, rmsLenKick, KickEntryHorizontal)
AddToField(Z, By, UNDULATOR_ZCENTER + 0.5 * DistanceBetweenKicks, rmsLenKick, KickExitHorizontal)
AddToField(Z, Bx, UNDULATOR_ZCENTER - 0.5 * DistanceBetweenKicks, rmsLenKick, KickEntryVertical)
AddToField(Z, Bx, UNDULATOR_ZCENTER + 0.5 * DistanceBetweenKicks, rmsLenKick, KickExitVertical)

# Field integrals
BxCorrI1 = IntegralVector(Z, Bx)
ByCorrI1 = IntegralVector(Z, By)
BxCorrI2 = IntegralVector(Z, BxCorrI1)
ByCorrI2 = IntegralVector(Z, ByCorrI1)
print "Before Correction 1st Integral Bx  %+.6E By %+.6E" % (FInt1Bx[-1],  FInt1By[-1])
print "After  Correction 1st Integral Bx  %+.6E By %+.6E" % (BxCorrI1[-1], ByCorrI1[-1])
print "Before Correction 2nd Integral Bx  %+.6E By %+.6E" % (FInt2Bx[-1],  FInt2By[-1])
print "After  Correction 2nd Integral Bx  %+.6E By %+.6E" % (BxCorrI2[-1], ByCorrI2[-1])







# Corrected field container
magFldCnt_Corr = SRWLMagFldC()
magFldCnt_Corr.allocate(1)

# Make magnetic field object for corrected field
magFldCnt_Corr.arMagFld[0] = SRWLMagFld3D( array('d', Bx), array('d', By), array('d', Bz), 1, 1, len(Z), 0.0, 0.0, Z[-1] - Z[0], 1, 1, None, None, _arZ=array('d', Z))

# Field interpolation method
magFldCnt_Corr.arMagFld[0].interp = 4

# ID center
magFldCnt_Corr.arXc[0] = 0.0
magFldCnt_Corr.arYc[0] = 0.0
magFldCnt_Corr.arZc[0] = UNDULATOR_ZCENTER


# Number of reps of field
magFldCnt_Corr.arMagFld[0].nRep = 1


# Get particle trajectory for corrected field
PartTraj_Corr = GetElectronTrajectory(magFldCnt_Corr, TRAJECTORY_START, TRAJECTORY_STOP)


# Get the Z Values for the plot
ZValues = numpy.linspace(PartTraj_Corr.ctStart, PartTraj_Corr.ctEnd, PartTraj_Corr.np)

# Graphs for trajectory from uncorrected data
gElectronX_Corr = TGraph( len(ZValues), array('d', ZValues), PartTraj_Corr.arX)
gElectronY_Corr = TGraph( len(ZValues), array('d', ZValues), PartTraj_Corr.arY)

gElectronX_Corr.SetTitle('Corrected Electron Trajectory in X')
gElectronX_Corr.GetXaxis().SetTitle('Z Position [m]')
gElectronX_Corr.GetYaxis().SetTitle('X Position [m]')
gElectronX_Corr.SetName("ElectronX_Corr")

gElectronY_Corr.SetTitle('Corrected Electron Trajectory in Y')
gElectronY_Corr.GetXaxis().SetTitle('Z Position [m]')
gElectronY_Corr.GetYaxis().SetTitle('Y Position [m]')
gElectronY_Corr.SetName("ElectronY_Corr")

gElectronX_Corr.Write()
gElectronY_Corr.Write()





# Define electron beam
ElecBeam = srwl_uti_src_e_beam('NSLS-II Low Beta Final')
ElecBeam.partStatMom1.x = 0.
ElecBeam.partStatMom1.y = 0.
ElecBeam.partStatMom1.z = TRAJECTORY_START
ElecBeam.partStatMom1.xp = 0
ElecBeam.partStatMom1.yp = 0



# Spectrum from Ideal field
[SpectrumIdealX, SpectrumIdealY] = GetUndulatorSpectrum(magFldCnt_Ideal, ElecBeam)
gSpectrumIdeal = TGraph( len(SpectrumIdealX), array('d', SpectrumIdealX), array('d', SpectrumIdealY) )
gSpectrumIdeal.SetName('Spectrum_Ideal')
gSpectrumIdeal.SetTitle('Simulated Spectrum')
gSpectrumIdeal.GetXaxis().SetTitle('Photon Energy [eV]')
gSpectrumIdeal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumIdeal.Write()


# Spectrum from corrected field
[SpectrumCorrX, SpectrumCorrY] = GetUndulatorSpectrum(magFldCnt_Corr, ElecBeam)
gSpectrumCorr = TGraph( len(SpectrumCorrX), array('d', SpectrumCorrX), array('d', SpectrumCorrY) )
gSpectrumCorr.SetName('Spectrum_DataCorr')
gSpectrumCorr.SetTitle('Simulated Spectrum from Field Measurements')
gSpectrumCorr.GetXaxis().SetTitle('Photon Energy [eV]')
gSpectrumCorr.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumCorr.Write()




StkP_Corr = SRWLStokes() # for power density
StkP_Corr.allocate(1, 100, 100) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
StkP_Corr.mesh.zStart = 30. #longitudinal position [m] at which power density has to be calculated
StkP_Corr.mesh.xStart = -0.02 #initial horizontal position [m]
StkP_Corr.mesh.xFin = 0.02 #final horizontal position [m]
StkP_Corr.mesh.yStart = -0.015 #initial vertical position [m]
StkP_Corr.mesh.yFin = 0.015 #final vertical position [m]
#StkP_Corr.mesh.eStart = 20000. #initial photon energy [eV]
#StkP_Corr.mesh.eFin = 20001. #final photon energy [eV]

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

print 'Begin power density calculation'
srwl.CalcPowDenSR(StkP_Corr, ElecBeam, 0, magFldCnt_Corr, arPrecP)

XValues = numpy.linspace(StkP_Corr.mesh.xStart, StkP_Corr.mesh.xFin, StkP_Corr.mesh.nx)
YValues = numpy.linspace(StkP_Corr.mesh.yStart, StkP_Corr.mesh.yFin, StkP_Corr.mesh.ny)
ZValues = numpy.array(StkP_Corr.arS[0:StkP_Corr.mesh.nx*StkP_Corr.mesh.ny]).reshape(StkP_Corr.mesh.ny, StkP_Corr.mesh.nx)
hPower_Corr = TH2F("PowerDensity_Corr", "Power Density", len(XValues), StkP_Corr.mesh.xStart, StkP_Corr.mesh.xFin, len(YValues), StkP_Corr.mesh.yStart, StkP_Corr.mesh.yFin)
for i in range( len(XValues) ):
  for j in range( len(YValues) ):
    hPower_Corr.SetBinContent(i+1, j+1, ZValues[i][j])
hPower_Corr.Write()




print 'Begin power density calculation Ideal'
srwl.CalcPowDenSR(StkP_Corr, ElecBeam, 0, magFldCnt_Ideal, arPrecP)

XValues = numpy.linspace(StkP_Corr.mesh.xStart, StkP_Corr.mesh.xFin, StkP_Corr.mesh.nx)
YValues = numpy.linspace(StkP_Corr.mesh.yStart, StkP_Corr.mesh.yFin, StkP_Corr.mesh.ny)
ZValues = numpy.array(StkP_Corr.arS[0:StkP_Corr.mesh.nx*StkP_Corr.mesh.ny]).reshape(StkP_Corr.mesh.ny, StkP_Corr.mesh.nx)
hPower_Ideal = TH2F("PowerDensity_Ideal", "Power Density", len(XValues), StkP_Corr.mesh.xStart, StkP_Corr.mesh.xFin, len(YValues), StkP_Corr.mesh.yStart, StkP_Corr.mesh.yFin)
for i in range( len(XValues) ):
  for j in range( len(YValues) ):
    hPower_Ideal.SetBinContent(i+1, j+1, ZValues[i][j])
hPower_Ideal.Write()



fROOT.Close()
