#!/usr/bin/env python
# Calculate the electron trajectory for measurements made in the lab

from SRWToolsUtil import *
import matplotlib.pyplot as plt

from uti_plot import *







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
plt.subplot(211)
plt.plot(ZValues, partTraj.arX)
plt.subplot(212)
plt.plot(ZValues, partTraj.arY)
plt.show()





