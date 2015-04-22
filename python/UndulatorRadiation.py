#!/usr/bin/env python
# Calculate the undulator ratiation

from SRWToolsUtil import *
import matplotlib.pyplot as plt





# grab the input date file from the command line
InFileName = sys.argv[1]



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
