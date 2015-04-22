# Path to SRW python lib and import as srw
import sys
sys.path.append('/home/dhidas/SRW/env/work/srw_python')
from srwlib import *



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




