#!/usr/bin/env python

import os
import numpy
from array import array
import sys

sys.argv.append('-b-')

from SRWToolsUtil import *


# extra root imports
from ROOT import TLine, TH2F










# Input and output files are given on command line
InFileName = sys.argv[1]

BaseFileName = os.path.splitext(os.path.basename(InFileName))[0]

# Open the input and output files.  Output file is ROOT format
fi = open(InFileName, 'r')
fROOT = TFile(BaseFileName + '.root', 'recreate')
fPEAKS_Ideal = open(BaseFileName + '.peaks.Ideal.dat', 'w')
#fPEAKS_Data = open(BaseFileName + '.peaks.Data.dat', 'w')
fPEAKS_Corr = open(BaseFileName + '.peaks.Corr.dat', 'w')

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

# Create arrays for magnetic field.  These contain the field data for the entire input file
Z  = array('d')
Bx = array('d')
By = array('d')
Bz = array('d')


# Loop over all entries in the input file
for l in fi:
  [tZ[0], tBx[0], tBy[0], tBz[0]] = map(float, l.split())

  # Make in SI unites and fill the Tree
  tZ[0] /= 1000.
  #t.Fill()

  # Fill program use arrays with data
  Z.append(tZ[0])
  Bx.append(tBx[0])
  By.append(tBy[0])
  Bz.append(tBz[0])

# Close file
fi.close()

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
gBy.SetMarkerStyle(3)
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
gMaxByChoppedAbs = TGraph(len(MaxByZChopped), array('d', MaxByZChopped), array('d', MaxByChoppedAbs))
gMaxByChoppedAbs.SetName('MaxByChoppedAbs')
gMaxByChoppedAbs.SetTitle('Calculated Max By')
gMaxByChoppedAbs.GetXaxis().SetTitle('Position [m]')
gMaxByChoppedAbs.GetYaxis().SetTitle('B_{Y} [T]')
fline = TF1('fline', 'pol1', -3000, 3000)
gMaxByChoppedAbs.Fit(fline, 'q')
gMaxByChoppedAbs.SetMarkerStyle(32)
gMaxByChoppedAbs.Write()
print 'Max By linear fit parameters [0]+[1]x::', fline.GetParameter(0), fline.GetParameter(1)


# Look at the pos vs neg peaks
cMax  = []
cMaxZ = []
cMin  = []
cMinZ = []
for i in range( len(MaxByChopped) ):
  x = MaxByChopped[i]
  if x >= 0:
    cMax.append(x)
    cMaxZ.append( MaxByZChopped[i] )
  else:
    cMin.append(x)
    cMinZ.append( MaxByZChopped[i] )









# Defining Magnetic Field container:
magFldCnt_Data = SRWLMagFldC()
magFldCnt_Data.allocate(1)


# Read data from file and make mag field object
magFldCnt_Data.arMagFld[0] = ReadHallProbeDataSRW(InFileName, UNDULATOR_ZSTART, UNDULATOR_ZEND)

# Field interpolation method
magFldCnt_Data.arMagFld[0].interp = 4

# ID center
magFldCnt_Data.arXc[0] = 0.0
magFldCnt_Data.arYc[0] = 0.0
magFldCnt_Data.arZc[0] = UNDULATOR_ZCENTER

# Number of reps of field
magFldCnt_Data.arMagFld[0].nRep = 1

# Get the electron trajectory
partTraj = GetElectronTrajectory(magFldCnt_Data, UNDULATOR_ZSTART, UNDULATOR_ZEND)

# Get the Z Values for the plot
ZValues = numpy.linspace(partTraj.ctStart, partTraj.ctEnd, partTraj.np)


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


# Define electron beam
elecBeam = srwl_uti_src_e_beam('NSLS-II Low Beta Final')
elecBeam.partStatMom1.x = 0.
elecBeam.partStatMom1.y = 0.
elecBeam.partStatMom1.z = -0.5 * PERIOD_LENGTH * (NPERIODS + 4)
elecBeam.partStatMom1.xp = 0
elecBeam.partStatMom1.yp = 0





# Get undulator and Ideal magnetic field based off of average peak By field
und = SRWLMagFldU([SRWLMagFldH(1, 'v', undBy, phBy, sBy, 1), SRWLMagFldH(1, 'h', undBx, phBx, sBx, 1)], undPer, numPer)
magFldCnt_Ideal = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements


stkF = SRWLStokes() #for spectral flux vs photon energy
stkF.allocate(10000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.zStart = 30. #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = 10. #initial photon energy [eV]
stkF.mesh.eFin = 60000. #final photon energy [eV]
stkF.mesh.xStart = -0.0015 #initial horizontal position [m]
stkF.mesh.xFin = 0.0015 #final horizontal position [m]
stkF.mesh.yStart = -0.00075 #initial vertical position [m]
stkF.mesh.yFin = 0.00075 #final vertical position [m]



#***********Precision Parameters
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = 21 #final UR harmonic to take into account
arPrecF[2] = 1.5 #longitudinal integration precision parameter
arPrecF[3] = 1.5 #azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

#*********************Calculation (SRWLIB function calls)
print('   Performing Spectral Flux (Stokes parameters) calculation ... ')
srwl.CalcStokesUR(stkF, elecBeam, und, arPrecF)
print('done')
ZValues = numpy.linspace(stkF.mesh.eStart, stkF.mesh.eFin, stkF.mesh.ne)
YValues = stkF.arS[0:stkF.mesh.ne]
gStokes = TGraph( len(ZValues), array('d', ZValues), array('d', YValues))
gStokes.SetName("Stokes")
gStokes.Write()




# Get the electron trajectory
partTraj_Ideal = GetElectronTrajectory(magFldCnt_Ideal, -1, 1)

# Get the Z Values for the plot
ZValues = numpy.linspace(partTraj_Ideal.ctStart, partTraj_Ideal.ctEnd, partTraj_Ideal.np)


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





[SpectrumIdealX, SpectrumIdealY] = GetUndulatorSpectrum(magFldCnt_Ideal, elecBeam)
gSpectrumIdeal = TGraph( len(SpectrumIdealX), array('d', SpectrumIdealX), array('d', SpectrumIdealY) )
gSpectrumIdeal.SetName('Spectrum_Ideal')
gSpectrumIdeal.SetTitle('Simulated Spectrum')
gSpectrumIdeal.GetXaxis().SetTitle('Photon Energy [eV]')
gSpectrumIdeal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumIdeal.Write()





# This is only here to test
#AddToField(Z, By, UNDULATOR_ZCENTER, 1.2, -0.00001)


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
sRangeX = (npAuxFldInt - 1) * StepSize(Z)

npAuxFldInt = len(By)
sRangeY = (npAuxFldInt - 1) * StepSize(Z)




# Distance between kicks
#DistanceBetweenKicks = PERIOD_LENGTH * (NPERIODS + 4)
DistanceBetweenKicks = UNDULATOR_ZEND - UNDULATOR_ZSTART
sCenX = 0
sCenY = 0

# RMS Kick Length in [m]
rmsLenKick = 0.050

KickEntryHorizontal =  0.5 * (sRangeY / DistanceBetweenKicks - 1) * auxI1Y - auxI2Y / DistanceBetweenKicks
KickExitHorizontal  = -0.5 * (sRangeY / DistanceBetweenKicks + 1) * auxI1Y + auxI2Y / DistanceBetweenKicks
KickEntryVertical   =  0.5 * (sRangeX / DistanceBetweenKicks - 1) * auxI1X - auxI2X / DistanceBetweenKicks
KickExitVertical    = -0.5 * (sRangeX / DistanceBetweenKicks + 1) * auxI1X + auxI2X / DistanceBetweenKicks


print 'DistanceBetweenKicks: ', DistanceBetweenKicks
print 'Kick EnH ExH EnV ExV: ', KickEntryHorizontal, KickExitHorizontal, KickEntryVertical, KickExitVertical


HalfDistBwKicks = 0.5 * DistanceBetweenKicks

AddToField(Z, By, UNDULATOR_ZCENTER - HalfDistBwKicks, rmsLenKick, KickEntryHorizontal)
AddToField(Z, By, UNDULATOR_ZCENTER + HalfDistBwKicks, rmsLenKick, KickExitHorizontal)
AddToField(Z, Bx, UNDULATOR_ZCENTER - HalfDistBwKicks, rmsLenKick, KickEntryVertical)
AddToField(Z, Bx, UNDULATOR_ZCENTER + HalfDistBwKicks, rmsLenKick, KickExitVertical)


#AddToField(Z, By, -sRangeY/2 + 3 * rmsLenKick, rmsLenKick, KickEntryHorizontal)
#AddToField(Z, By,  sRangeY/2 - 3 * rmsLenKick, rmsLenKick, KickExitHorizontal)
#AddToField(Z, Bx, -sRangeX/2 + 3 * rmsLenKick, rmsLenKick, KickEntryVertical)
#AddToField(Z, Bx,  sRangeX/2 - 3 * rmsLenKick, rmsLenKick, KickExitVertical)


BxCorrI1 = IntegralVector(Z, Bx)
ByCorrI1 = IntegralVector(Z, By)
BxCorrI2 = IntegralVector(Z, BxCorrI1)
ByCorrI2 = IntegralVector(Z, ByCorrI1)
print "Befor Correction 1st Integral Bx By", FInt1Bx[-1],  FInt1By[-1]
print "After Correction 1st Integral Bx By", BxCorrI1[-1], ByCorrI1[-1]
print "Befor Correction 2nd Integral Bx By", FInt2Bx[-1],  FInt2By[-1]
print "After Correction 2nd Integral Bx By", BxCorrI2[-1], ByCorrI2[-1]




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


# Get the electron trajectory
#partTraj_Corr = GetElectronTrajectory(magFldCnt_Corr, UNDULATOR_ZSTART - 1.50, UNDULATOR_ZEND + 1.50)
partTraj_Corr = GetElectronTrajectory(magFldCnt_Corr, -1, 1)

# Get the Z Values for the plot
ZValues_Corr = numpy.linspace(partTraj_Corr.ctStart, partTraj_Corr.ctEnd, partTraj_Corr.np)


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



#[SpectrumDataX, SpectrumDataY] = GetUndulatorSpectrum(magFldCnt_Data, elecBeam)
#gSpectrumData = TGraph( len(SpectrumDataX), array('d', SpectrumDataX), array('d', SpectrumDataY) )
#gSpectrumData.SetName('Spectrum_Data')
#gSpectrumData.SetTitle('Simulated Spectrum from Field Measurements')
#gSpectrumData.GetXaxis().SetTitle('Photon Energy [eV]')
#gSpectrumData.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
#gSpectrumData.Write()

stkP = SRWLStokes() #for power density
stkP.allocate(1, 100, 100) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = 30. #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -0.02 #initial horizontal position [m]
stkP.mesh.xFin = 0.02 #final horizontal position [m]
stkP.mesh.yStart = -0.015 #initial vertical position [m]
stkP.mesh.yFin = 0.015 #final vertical position [m]
#stkP.mesh.eStart = 20000. #initial photon energy [eV]
#stkP.mesh.eFin = 20001. #final photon energy [eV]

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

print 'Begin power density calculation'
srwl.CalcPowDenSR(stkP, elecBeam, 0, magFldCnt_Corr, arPrecP)

XValues = numpy.linspace(stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx)
YValues = numpy.linspace(stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny)
ZValues = numpy.array(stkP.arS[0:stkP.mesh.nx*stkP.mesh.ny]).reshape(stkP.mesh.ny, stkP.mesh.nx)
hPower = TH2F("PowerDensity", "Power Density", len(XValues), stkP.mesh.xStart, stkP.mesh.xFin, len(YValues), stkP.mesh.yStart, stkP.mesh.yFin)
for i in range( len(XValues) ):
  for j in range( len(YValues) ):
    hPower.SetBinContent(i+1, j+1, ZValues[i][j])
hPower.Write()
#gStokes.SetName("Stokes")
#gStokes.Write()


exit(0)

[SpectrumCorrX, SpectrumCorrY] = GetUndulatorSpectrum(magFldCnt_Corr, elecBeam)
gSpectrumCorr = TGraph( len(SpectrumCorrX), array('d', SpectrumCorrX), array('d', SpectrumCorrY) )
gSpectrumCorr.SetName('Spectrum_DataCorr')
gSpectrumCorr.SetTitle('Simulated Spectrum from Field Measurements')
gSpectrumCorr.GetXaxis().SetTitle('Photon Energy [eV]')
gSpectrumCorr.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
gSpectrumCorr.Write()




# Run peak finding on all spectra
hSpectrumIdeal = TGraphToTH1F(gSpectrumIdeal)
#hSpectrumData = TGraphToTH1F(gSpectrumData)
hSpectrumCorr = TGraphToTH1F(gSpectrumCorr)
hSpectrumIdeal.Write()
#hSpectrumData.Write()
hSpectrumCorr.Write()
PeaksIdeal = FindPeaksInHistogram(hSpectrumIdeal)
#PeaksData = FindPeaksInHistogram(hSpectrumData)
PeaksCorr = FindPeaksInHistogram(hSpectrumCorr)

# Show the peaks on the histogram plot
GetCanvasWithHistAndPeakMarkers(hSpectrumIdeal, PeaksIdeal, BaseFileName + '_Spectrum_Ideal' ).Write()
GetCanvasWithHistAndPeakMarkers(hSpectrumCorr, PeaksCorr,BaseFileName + '_Spectrum_Corr' ).Write()
#GetCanvasWithHistAndPeakMarkers(hSpectrumData, PeaksData,BaseFileName + '_Spectrum_Data' ).Write()


PeaksIdealSorted = []
#PeaksDataSorted = []
PeaksCorrSorted = []
for peak in PeaksIdeal:
  PeaksIdealSorted.append(peak)
#for peak in PeaksData:
#  PeaksDataSorted.append(peak)
for peak in PeaksCorr:
  PeaksCorrSorted.append(peak)
for peak in PeaksIdeal:
  PeaksIdealSorted.append(peak)

# Now do sorting for all three
PeaksCorrSorted.sort()
#PeaksDataSorted.sort()
PeaksIdealSorted.sort()

for i in range( len(PeaksCorrSorted) ):
  p = PeaksCorrSorted[i]
  fPEAKS_Corr.write('%2i  %8.1f  %.6E\n' % (i, p, PeaksCorr[p]))
#for i in range( len(PeaksDataSorted) ):
#  p = PeaksDataSorted[i]
#  fPEAKS_Data.write('%2i %8.1f  %.6E\n' % (i, p, PeaksData[p]))
for i in range( len(PeaksIdealSorted) ):
  p = PeaksIdealSorted[i]
  fPEAKS_Ideal.write('%2i  %8.1f  %.6E\n' % (i, p, PeaksIdeal[p]))



fPEAKS_Corr.close()
#fPEAKS_Data.close()
fPEAKS_Ideal.close()











fROOT.Write()
fROOT.Close()

