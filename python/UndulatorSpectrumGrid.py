#!/usr/bin/env python

import os
import numpy
from array import array
import sys

sys.argv.append('-b-')

from SRWToolsUtil import *


# extra root imports
from ROOT import TLine, TH1F



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

BaseFileName = os.path.splitext(os.path.basename(InFileName))[0]

# Name of output file for this section
OutFileNameData = BaseFileName + '_Section_' + str(SectionNumber) + '_Spectrum.dat'
print 'Output file for spectrum calculation: ', OutFileNameData
OutFileData = open(OutFileNameData, 'w')
OutFileData.write(str(RandomSeed) + '\n')

# Open the input and output files.  Output file is ROOT format
fi = open(InFileName, 'r')
fROOT = None
if (DEBUG):
  fROOT = TFile(BaseFileName + '_Section_' + str(SectionNumber) + '_Spectrum.root', 'recreate')
  fROOT.cd()

# Create variables for TTree filling
tZ  = float()
tBx = float()
tBy = float()
tBz = float()


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
[MaxByZ, MaxBy] = FindMinMaxFromFit(MaxListInd, Z, By, fROOT)


if (DEBUG):
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


# Peak By to use in calculations
PEAK_BY = AvgMaxBy



if (DEBUG):
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

if (DEBUG):
  gMaxByChoppedAbs.Write()


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



elecX0 = elecBeam.partStatMom1.x
elecXp0 = elecBeam.partStatMom1.xp
elecY0 = elecBeam.partStatMom1.y
elecYp0 = elecBeam.partStatMom1.yp
elecGamma0 = elecBeam.partStatMom1.gamma
elecE0 = elecGamma0*(0.51099890221e-03) #Assuming electrons 

elecSigXe2 = elecBeam.arStatMom2[0] #<(x-x0)^2>
elecMXXp = elecBeam.arStatMom2[1] #<(x-x0)*(xp-xp0)>
elecSigXpe2 = elecBeam.arStatMom2[2] #<(xp-xp0)^2>
elecSigYe2 =elecBeam.arStatMom2[3] #<(y-y0)^2>
elecMYYp = elecBeam.arStatMom2[4] #<(y-y0)*(yp-yp0)>
elecSigYpe2 = elecBeam.arStatMom2[5] #<(yp-yp0)^2>
elecRelEnSpr = sqrt(elecBeam.arStatMom2[10]) #<(E-E0)^2>/E0^2
elecAbsEnSpr = elecE0*elecRelEnSpr
#print('DEBUG MESSAGE: elecAbsEnSpr=', elecAbsEnSpr)

multX = 0.5/(elecSigXe2*elecSigXpe2 - elecMXXp*elecMXXp)
BX = elecSigXe2*multX
GX = elecSigXpe2*multX
AX = elecMXXp*multX
SigPX = 1/sqrt(2*GX)
SigQX = sqrt(GX/(2*(BX*GX - AX*AX)))
multY = 0.5/(elecSigYe2*elecSigYpe2 - elecMYYp*elecMYYp)
BY = elecSigYe2*multY
GY = elecSigYpe2*multY
AY = elecMYYp*multY
SigPY = 1/sqrt(2*GY)
SigQY = sqrt(GY/(2*(BY*GY - AY*AY)))


# Get undulator and Ideal magnetic field based off of average peak By field
und = SRWLMagFldU([SRWLMagFldH(1, 'v', undBy, phBy, sBy, 1), SRWLMagFldH(1, 'h', undBx, phBx, sBx, 1)], undPer, numPer)
magFldCnt_Ideal = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements


# Histogram 
hX  = TH1F('partStatMom1.x',  'partStatMom1.x',  100, -1e-4, 1e-4)
hY  = TH1F('partStatMom1.y',  'partStatMom1.y',  100, -1e-5, 1e-5)
hXP = TH1F('partStatMom1.xp', 'partStatMom1.xp', 100, -1e-4, 1e-4)
hYP = TH1F('partStatMom1.yp', 'partStatMom1.yp', 100, -1e-5, 1e-5)
hG = TH1F('artStatMom1.gamma', 'artStatMom1.gamma', 100, 5700, 6100)

SpectrumAverages_Ideal = []
SpectrumXValues_Ideal = []

for i in range(5):
  print 'Electron number in this section:', i

  elecBeamCopy = deepcopy(elecBeam)

  auxPXp = SigQX * random.gauss(0, 1)
  auxPX = SigPX * random.gauss(0, 1) + AX * auxPXp / GX
  elecBeamCopy.partStatMom1.x = elecX0 + auxPX
  elecBeamCopy.partStatMom1.xp = elecXp0 + auxPXp
  auxPYp = SigQY * random.gauss(0, 1)
  auxPY = SigPY * random.gauss(0, 1) + AY * auxPYp / GY
  elecBeamCopy.partStatMom1.y = elecY0 + auxPY
  elecBeamCopy.partStatMom1.yp = elecYp0 + auxPYp
  elecBeamCopy.partStatMom1.gamma = elecGamma0 * (1 + elecAbsEnSpr * random.gauss(0, 1) / elecE0)

  if (DEBUG):
    hX.Fill(elecBeamCopy.partStatMom1.x)
    hY.Fill(elecBeamCopy.partStatMom1.y)
    hXP.Fill(elecBeamCopy.partStatMom1.xp)
    hYP.Fill(elecBeamCopy.partStatMom1.yp)
    hG.Fill(elecBeamCopy.partStatMom1.gamma)

  print elecBeamCopy.partStatMom1.x

  # Get the spectrum
  if SectionNumber == 0:
    [X, Y] = GetUndulatorSpectrum(magFldCnt_Ideal, elecBeam)
  else:
    [X, Y] = GetUndulatorSpectrum(magFldCnt_Ideal, elecBeamCopy)

  if i == 0:
    SpectrumXValues_Ideal = X[:]

  AddToRunningAverages (SpectrumAverages_Ideal, i, Y)



  if (DEBUG):
    gSpectrumIdeal = TGraph( len(X), array('d', X), array('d', Y) )
    gSpectrumIdeal.SetName('Spectrum_Ideal_' + str(i))
    gSpectrumIdeal.SetTitle('Simulated Spectrum')
    gSpectrumIdeal.GetXaxis().SetTitle('Photon Energy [eV]')
    gSpectrumIdeal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
    gSpectrumIdeal.Write()

  if SectionNumber == 0:
    break



if (DEBUG):
  gSpectrumIdeal = TGraph( len(SpectrumXValues_Ideal), array('d', SpectrumXValues_Ideal), array('d', SpectrumAverages_Ideal) )
  gSpectrumIdeal.SetName('Spectrum_Ideal_Average')
  gSpectrumIdeal.SetTitle('Simulated Spectrum')
  gSpectrumIdeal.GetXaxis().SetTitle('Photon Energy [eV]')
  gSpectrumIdeal.GetYaxis().SetTitle('Intensity photons/s/.1%bw/mm^{2}')
  gSpectrumIdeal.Write()
  hX.Write()
  hY.Write()
  hXP.Write()
  hYP.Write()
  hG.Write()
  







for i in range( len(SpectrumXValues_Ideal) ):
  OutFileData.write( str(SpectrumXValues_Ideal[i]) + ' ' + str(SpectrumAverages_Ideal[i]) + '\n' )



exit(0)



if (DEBUG):
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








if (DEBUG):
  fROOT.Close()





print 'done.'
