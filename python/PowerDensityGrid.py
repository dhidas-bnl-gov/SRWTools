#!/usr/bin/env python

import os
import numpy
from array import array
import sys

sys.argv.append('-b-')

from SRWToolsUtil import *


# extra root imports
from ROOT import TLine, TH1F, TH2F



# For debug output and more plots
DEBUG = True







# Input and output files are given on command line
InFileName = sys.argv[1]
SectionNumber = int(sys.argv[2])

print 'Input file:', InFileName
print 'Section number:', SectionNumber

# Get a random seed based on section number
RandomSeed = 131147 + SectionNumber
random.seed(RandomSeed)
print 'Set random.seed() to', RandomSeed


# Name of output files for this section
BaseFileName = os.path.splitext(os.path.basename(InFileName))[0]
OutFileNameDataIdeal = BaseFileName + '_PowerDensityIdeal_Section_' + str(SectionNumber).zfill(4) + '.dat' if SectionNumber > 0 else BaseFileName + '_PowerDensityIdeal_SingleElectron.dat'
OutFileNameDataCorr = BaseFileName + '_PowerDensityCorr_Section_' + str(SectionNumber).zfill(4) + '.dat' if SectionNumber > 0 else BaseFileName + '_PowerDensityCorr_SingleElectron.dat'
OutFileNameRoot = BaseFileName + '_PowerDensity_Section_' + str(SectionNumber).zfill(4) + '.root' if SectionNumber > 0 else BaseFileName + '_PowerDensity_SingleElectron.root'
print 'Output file for spectrum calculation: ', OutFileNameDataIdeal, OutFileNameDataCorr, OutFileNameRoot

# Open output file and print seed as the first line (seed checking is done later in combination)
fOutFileDataIdeal = open(OutFileNameDataIdeal, 'w')
fOutFileDataCorr = open(OutFileNameDataCorr, 'w')
fOutFileDataIdeal.write(str(RandomSeed) + '\n')
fOutFileDataCorr.write(str(RandomSeed) + '\n')

# Open the input and output ROOT file.
fi = open(InFileName, 'r')
fROOT = TFile(OutFileNameRoot, 'recreate')
fROOT.cd()


# Create arrays for magnetic field.  These contain the field data for the entire input file
Z  = array('d')
Bx = array('d')
By = array('d')
Bz = array('d')


# Loop over all entries in the input file
for l in fi:
  [tZ, tBx, tBy, tBz] = map(float, l.split())

  # Make in SI unites and fill the Tree
  tZ /= 1000.

  # Fill program use arrays with data
  Z.append(tZ)
  Bx.append(tBx)
  By.append(tBy)
  Bz.append(tBz)

# Close file
fi.close()

# Get the max and mins in the field
[MaxListInd, MaxListBy] = FindMaxAndMins(Z, By)
[MaxByZ, MaxBy] = FindMinMaxFromFit(MaxListInd, Z, By, None)





if SectionNumber == 0:
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


# Draw Undulator fields and START, CENTER, END
if SectionNumber == 0:
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

if SectionNumber == 0:
  gMaxByChoppedAbs.Write()





# Field integrals
FInt1Bx = IntegralVector(Z, Bx)
FInt2Bx = IntegralVector(Z, FInt1Bx)
FInt1By = IntegralVector(Z, By)
FInt2By = IntegralVector(Z, FInt1By)

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


print 'DistanceBetweenKicks: ', DistanceBetweenKicks
print 'Kick EnH ExH EnV ExV: ', KickEntryHorizontal, KickExitHorizontal, KickEntryVertical, KickExitVertical


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



elecX0     = elecBeam.partStatMom1.x
elecXp0    = elecBeam.partStatMom1.xp
elecY0     = elecBeam.partStatMom1.y
elecYp0    = elecBeam.partStatMom1.yp
elecGamma0 = elecBeam.partStatMom1.gamma
elecE0     = elecGamma0*(0.51099890221e-03) #Assuming electrons 

elecSigXe2   = elecBeam.arStatMom2[0] #<(x-x0)^2>
elecMXXp     = elecBeam.arStatMom2[1] #<(x-x0)*(xp-xp0)>
elecSigXpe2  = elecBeam.arStatMom2[2] #<(xp-xp0)^2>
elecSigYe2   = elecBeam.arStatMom2[3] #<(y-y0)^2>
elecMYYp     = elecBeam.arStatMom2[4] #<(y-y0)*(yp-yp0)>
elecSigYpe2  = elecBeam.arStatMom2[5] #<(yp-yp0)^2>
elecRelEnSpr = sqrt(elecBeam.arStatMom2[10]) #<(E-E0)^2>/E0^2
elecAbsEnSpr = elecE0*elecRelEnSpr

multX = 0.5 / (elecSigXe2 * elecSigXpe2 - elecMXXp * elecMXXp)
BX    = elecSigXe2 * multX
GX    = elecSigXpe2 * multX
AX    = elecMXXp * multX
SigPX = 1 / sqrt(2 * GX)
SigQX = sqrt(GX / (2 * (BX * GX - AX * AX)))
multY = 0.5 / (elecSigYe2 * elecSigYpe2 - elecMYYp * elecMYYp)
BY    = elecSigYe2 * multY
GY    = elecSigYpe2 * multY
AY    = elecMYYp * multY
SigPY = 1 / sqrt(2 * GY)
SigQY = sqrt(GY / (2 * (BY * GY - AY * AY)))


# Get undulator and Ideal magnetic field based off of average peak By field
und = SRWLMagFldU([SRWLMagFldH(1, 'v', undBy, phBy, sBy, 1), SRWLMagFldH(1, 'h', undBx, phBx, sBx, 1)], undPer, numPer)
magFldCnt_Ideal = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements


# Histogram 
hX  = TH1F('partStatMom1.x',      'partStatMom1.x',      100, -1e-4, 1e-4)
hY  = TH1F('partStatMom1.y',      'partStatMom1.y',      100, -1e-5, 1e-5)
hXP = TH1F('partStatMom1.xp',     'partStatMom1.xp',     100, -1e-4, 1e-4)
hYP = TH1F('partStatMom1.yp',     'partStatMom1.yp',     100, -1e-5, 1e-5)
hG  = TH1F('partStatMom1.energy', 'partStatMom1.energy', 100,   2.9, 3.1)



# Variables for power density calculation
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

# Grid values in X and Y
XValues = numpy.linspace(stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx)
YValues = numpy.linspace(stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny)

PowerDensityXY_Ideal = []
PowerDensityXY_Corr  = []


for i in range(3):
  print 'Electron number in this section:', i

  # Copy the electron beam for smearing
  elecBeamCopy = deepcopy(elecBeam)

  # Smear X and X', Y, Y', and gamma
  auxPXp = SigQX * random.gauss(0, 1)
  auxPX  = SigPX * random.gauss(0, 1) + AX * auxPXp / GX
  elecBeamCopy.partStatMom1.x  = elecX0 + auxPX
  elecBeamCopy.partStatMom1.xp = elecXp0 + auxPXp
  auxPYp = SigQY * random.gauss(0, 1)
  auxPY  = SigPY * random.gauss(0, 1) + AY * auxPYp / GY
  elecBeamCopy.partStatMom1.y  = elecY0 + auxPY
  elecBeamCopy.partStatMom1.yp = elecYp0 + auxPYp
  elecBeamCopy.partStatMom1.gamma = elecGamma0 * (1 + elecAbsEnSpr * random.gauss(0, 1) / elecE0)

  # Fill histograms for diagnostics
  hX.Fill(elecBeamCopy.partStatMom1.x)
  hY.Fill(elecBeamCopy.partStatMom1.y)
  hXP.Fill(elecBeamCopy.partStatMom1.xp)
  hYP.Fill(elecBeamCopy.partStatMom1.yp)
  hG.Fill(elecBeamCopy.partStatMom1.gamma * (0.51099890221e-03))


  # Get the spectrum for Corrected field
  if SectionNumber == 0:
    srwl.CalcPowDenSR(stkP, elecBeam, 0, magFldCnt_Ideal, arPrecP)
  else:
    srwl.CalcPowDenSR(stkP, elecBeamCopy, 0, magFldCnt_Ideal, arPrecP)
  ZValues = numpy.array(stkP.arS[0:stkP.mesh.nx*stkP.mesh.ny]).reshape(stkP.mesh.ny, stkP.mesh.nx)
  PowerDensityXY_Ideal = AddToRunningAverages2D(PowerDensityXY_Ideal, i, ZValues)


  # Get the spectrum for Corrected field
  if SectionNumber == 0:
    srwl.CalcPowDenSR(stkP, elecBeam, 0, magFldCnt_Corr, arPrecP)
  else:
    srwl.CalcPowDenSR(stkP, elecBeamCopy, 0, magFldCnt_Corr, arPrecP)
  ZValues = numpy.array(stkP.arS[0:stkP.mesh.nx*stkP.mesh.ny]).reshape(stkP.mesh.ny, stkP.mesh.nx)
  PowerDensityXY_Corr = AddToRunningAverages2D(PowerDensityXY_Corr, i, ZValues)


  if SectionNumber == 0:
    break





# need to implement X, Y limits
if DEBUG or SectionNumber == 0:
  hPowerDensity_Ideal = TH2F("PowerDensity_Ideal", "Power Density Ideal Field", len(PowerDensityXY_Ideal), stkP.mesh.xStart, stkP.mesh.xFin, len(PowerDensityXY_Ideal[0]), stkP.mesh.yStart, stkP.mesh.yFin)
  for j in range( len(PowerDensityXY_Ideal) ):
    for k in range( len(PowerDensityXY_Ideal[j]) ):
      hPowerDensity_Ideal.SetBinContent(j+1, k+1, PowerDensityXY_Ideal[j][k])
  hPowerDensity_Ideal.SetTitle('Simulated Spectrum')
  hPowerDensity_Ideal.GetXaxis().SetTitle('Photon Energy [eV]')
  hPowerDensity_Ideal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
  hPowerDensity_Ideal.Write()

  hPowerDensity_Corr = TH2F("PowerDensity_Corr", "Power Density Corrected Field", len(PowerDensityXY_Corr), stkP.mesh.xStart, stkP.mesh.xFin, len(PowerDensityXY_Corr[0]), stkP.mesh.yStart, stkP.mesh.yFin)
  for j in range( len(PowerDensityXY_Corr) ):
    for k in range( len(PowerDensityXY_Corr[j]) ):
      hPowerDensity_Corr.SetBinContent(j+1, k+1, PowerDensityXY_Corr[j][k])
  hPowerDensity_Corr.SetTitle('Simulated Spectrum')
  hPowerDensity_Corr.GetXaxis().SetTitle('Photon Energy [eV]')
  hPowerDensity_Corr.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
  hPowerDensity_Corr.Write()




for j in range( len(PowerDensityXY_Ideal) ):
  for k in range( len(PowerDensityXY_Ideal[j]) ):
    fOutFileDataIdeal.write('%+.9E ' % PowerDensityXY_Ideal[j][k])
  fOutFileDataIdeal.write('\n')

for j in range( len(PowerDensityXY_Corr) ):
  for k in range( len(PowerDensityXY_Corr[j]) ):
    fOutFileDataCorr.write('%+.9E ' % PowerDensityXY_Corr[j][k])
  fOutFileDataCorr.write('\n')

fOutFileDataIdeal.close()
fOutFileDataCorr.close()




if (DEBUG and SectionNumber != 0):
  hX.Write()
  hY.Write()
  hXP.Write()
  hYP.Write()
  hG.Write()
  

# Close root file
fROOT.Close()




print 'done.'
