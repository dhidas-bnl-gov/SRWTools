#!/usr/bin/env python
# Calculate the electron trajectory and undulator radiation emitted for measurements
# made in the lab using SRW

from SRWToolsUtil import *

import sys
import numpy as np
import matplotlib.pyplot as plt






# grab the input date file from the command line
InFileName = sys.argv[1]


# this is the particle which will travel through the field
part = SRWLParticle()
part.x     =  0.0 # Initial position x [m]
part.y     =  0.0 # Initial position y [m]
part.xp    =  0.0 # Iitial velocity
part.yp    =  0.0 # Initial velocity
part.gamma =  3/0.51099890221e-03 # Relative Energy
part.relE0 =  1 # Electron Rest Mass
part.nq    = -1 # Electron Charge



# Defining Magnetic Field container:
magFldCnt = SRWLMagFldC()
magFldCnt.allocate(1)


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

# Starting position in z of particle
part.z = 0.0 - 0.5*magFldCnt.arMagFld[0].rz

# Trajectory structure, where the results will be stored
partTraj = SRWLPrtTrj()
partTraj.partInitCond = part

# Number of trajectory points
partTraj.allocate(10001, True)
partTraj.ctStart = 0
partTraj.ctEnd = magFldCnt.arMagFld[0].rz


# Calculation (SRWLIB function call)
partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, [1])

# Get the Z Values for the plot
ZValues = [float(x) * ((partTraj.ctEnd - partTraj.ctStart) / float(partTraj.np)) for x in range(0, partTraj.np)]

# Convert to mm from m
for i in range(partTraj.np):
    partTraj.arX[i] *= 1000
    partTraj.arY[i] *= 1000
 
 
# Draw some plots
#plt.subplot(211)
#plt.plot(ZValues, partTraj.arX)
#plt.subplot(212)
#plt.plot(ZValues, partTraj.arY)
#plt.show()





#***********Electron Beam
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



#***********Precision
meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.01 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]

#**********************Calculation (SRWLIB function calls)
srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)


print [wfr1.mesh.eStart, wfr1.mesh.eFin, wfr1.mesh.ne]
XValues = [(wfr1.mesh.eStart + float(x) * (wfr1.mesh.eFin - wfr1.mesh.eStart) / wfr1.mesh.ne) for x in range(wfr1.mesh.ne)]
plt.plot(XValues, arI1)
plt.show()
