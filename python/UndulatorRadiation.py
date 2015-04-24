#!/usr/bin/env python
# Calculate the undulator ratiation

from SRWToolsUtil import *
import matplotlib.pyplot as plt






# grab the input date file from the command line
InFileName = ''
if len(sys.argv) == 2:
  InFileName = sys.argv[1]
elif len(sys.argv) > 2:
  print 'Usage: ', sys.argv[0], '[InFileName]'
  print 'If InFileName not specified, will run for simulated undulator'
  exit(0)



# Defining Magnetic Field container:
magFldCnt = SRWLMagFldC()
magFldCnt.allocate(1)



if InFileName == '':
  # define the undulator specs
  numPer = 70
  undPer = 0.021
  Bx = 0.0
  By = 1.2
  phBx = 0
  phBy = 0
  sBx = 1
  sBy = 1
  xcID = 0
  ycID = 0
  zcID = 0

  und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1), SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer, numPer)
  magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements
else :
  # Read data from file and make mag field object
  magFldCnt.arMagFld[0] = ReadHallProbeDataSRW(InFileName)

  # Field interpolation method
  magFldCnt.arMagFld[0].interp = 4

  # ID center
  magFldCnt.arXc[0] = 0.0
  magFldCnt.arYc[0] = 0.0


  # Number of reps of field
  magFldCnt.arMagFld[0].nRep = 1

  # Center in Z of ID
  magFldCnt.arZc[0] = 0.0



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

# Make some plots
plt.plot(XValues, arI1)
plt.show()
