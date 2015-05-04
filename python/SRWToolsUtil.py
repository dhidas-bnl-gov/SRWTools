# Path to SRW python lib and import as srw
import sys
sys.path.append('/home/dhidas/SRW/env/work/srw_python')
sys.path.append('/Users/dhidas/SRW/env/work/srw_python')
from srwlib import *

from ROOT import TFile, TTree, TGraph, TCanvas, TF1


def ReadHallProbeData (InFileName, ZMin = -999, ZMax = 999):
  "Read a data file from the teststand and return a list of 3 lists of floats.  The convention is [mm] and [T]"

  f = open(InFileName, 'r')

  Z  = []
  Bx = []
  By = []
  Bz = []

  for l in f:
    [z, bx, by, bz] = map(float, l.split())
    z /= 1000.
    if z > ZMin and z < ZMax:
      Z.append(z)
      Bx.append(bx)
      By.append(by)
      Bz.append(bz)
    

  return [Z, Bx, By, Bz]






def ReadHallProbeDataSRW (InFileName, ZMin = -0.850, ZMax = 0.850):
  "Read a data file from the teststand and return list of parameters needed for simulation"

  f = open(InFileName, 'r')

  Z  = []
  Bx = []
  By = []
  Bz = []

  npXY = 1
  npZ  = 0

  for l in f:
    [z, bx, by, bz] = map(float, l.split())
    z /= 1000.
    if z > ZMin and z < ZMax:
      Z.append(z)
      Bx.append(bx)
      By.append(by)
      Bz.append(bz)
      npZ += 1


  ZStep = (ZMax - ZMin) / npZ

  locArZ  = array('d', [0]*npZ)
  locArBx = array('d', [0]*npZ)
  locArBy = array('d', [0]*npZ)
  locArBz = array('d', [0]*npZ)

  for i in range(npZ):
    locArZ[i]  =  Z[i]
    locArBx[i] = Bx[i]
    locArBy[i] = By[i]
    locArBz[i] = Bz[i]
    

  return SRWLMagFld3D(locArBx, locArBy, locArBz, npXY, npXY, npZ, 0.0, 0.0, (npZ)*ZStep, 1, 1, None, None, _arZ=locArZ)






def IntegralVector (X, Y):
  "Integrate using simple trapazoid rule.  return vector: integral up to each point"

  if len(X) != len(Y):
    print 'ERROR: IntegralVector (X, Y) lengths not the same'
     
     

  I = [0.0]
  for i in range(1, len(X)):
    dx = X[i] - X[i-1]
    trap = dx * (Y[i-1] + 0.5 * (Y[i] - Y[i-1]))
    I.append(I[-1:][0] + trap)

  return I







def GetUndulatorSpectrum (magFldCnt,):
  "Get the spectrum given und"

  # Electron Beam
  elecBeam = srwl_uti_src_e_beam('NSLS-II Low Beta Final')
  elecBeam.partStatMom1.x = 0.
  elecBeam.partStatMom1.y = 0.
  elecBeam.partStatMom1.z = -1.8 #-0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
  elecBeam.partStatMom1.xp = 0
  elecBeam.partStatMom1.yp = 0

  # For spectrum vs photon energy
  wfr1 = SRWLWfr()
  wfr1.allocate(10000, 1, 1)
  wfr1.mesh.zStart = 20.
  wfr1.mesh.eStart = 10.
  wfr1.mesh.eFin = 30000.
  wfr1.mesh.xStart = 0.
  wfr1.mesh.xFin = 0
  wfr1.mesh.yStart = 0
  wfr1.mesh.yFin = 0
  wfr1.partBeam = elecBeam



  # Precision
  meth = 1
  relPrec = 0.01
  zStartInteg = 0
  zEndInteg = 0
  npTraj = 20000
  useTermin = 1
  sampFactNxNyForProp = 0
  arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]

  # Calculation (SRWLIB function calls)
  srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
  arI1 = array('f', [0]*wfr1.mesh.ne)
  srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)

  # Get X-values for plotting
  XValues = [(wfr1.mesh.eStart + float(x) * (wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne) for x in range(wfr1.mesh.ne)]

  return [XValues, arI1]





def FindMaxAndMins (Z, By):
  "find the max and mins in the By distribution"

  # Loop over data and find min/max for every point
  IsAbove = False
  IsBelow = False
  MaxBy = 0
  MinBy = 0
  ZMax  = 0
  ZMin  = 0
  ZMaxI = 0
  ZMinI = 0
  MaxListInd  = []
  MaxListBy   = []

  BThreshold = 0.0004
  for i in range( len(Z) ):
    if By[i] < -BThreshold:
      if IsAbove:
        MaxListInd.append(ZMaxI)
        MaxListBy.append(MaxBy)
        MaxBy = 0
      IsBelow = True
      IsAbove = False
    if By[i] > BThreshold:
      if IsBelow:
        MaxListInd.append(ZMinI)
        MaxListBy.append(MinBy)
        MinBy = 0
      IsBelow = False
      IsAbove = True

    if IsAbove and By[i] > MaxBy:
      MaxBy = By[i]
      ZMax  = Z[i]
      ZMaxI = i
    if IsBelow and By[i] < MinBy:
      MinBy = By[i]
      ZMin  = Z[i]
      ZMinI = i
  if IsAbove:
    MaxListInd.append(ZMaxI)
    MaxListBy.append(MaxBy)
  if IsBelow:
    MaxListInd.append(ZMinI)
    MaxListBy.append(MinBy)

  print 'Number of max/min seen I: ', len(MaxListInd)

  return [MaxListInd, MaxListBy]




def FindMinMaxFromFit (MaxListInd, Z, By):
  "Return the max B and Z locations from a fit"

  # Where is the calculated max By from the fit
  MaxBy  = []
  MaxByZ = []

  # Calculate the max based on pol2 fit maximum.  Assumption is to use max +-NFitWidth
  NFitWidth = 1
  for i in MaxListInd:
    x = []
    y = []

    for j in range(-NFitWidth, NFitWidth + 1):
      x.append(Z[i+j])
      if By[i] >= 0:
        y.append(By[i+j])
      else:
        y.append(-By[i+j])

    g = TGraph(NFitWidth*2 + 1, array('d', x), array('d', y))
    g.SetName('gFit' + str(i))
    FitFunction = TF1("FitFunction_"+str(i), "pol2", -4, 4)
    g.Fit(FitFunction, 'q')
    g.Write()
    MaxByZ.append(FitFunction.GetMaximumX())
    if By[i] >= 0:
      MaxBy.append(FitFunction.GetMaximum())
    else:
      MaxBy.append(-FitFunction.GetMaximum())

  return [MaxByZ, MaxBy]
