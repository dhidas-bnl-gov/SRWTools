#from srwlib import *



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
      npZ += 1



    

  return [Z, Bx, By, Bz]






def ReadHallProbeDataSRW (InFileName, ZMin = -999, ZMax = 999):
  "Read a data file from the teststand and return SRW MagFld3D object"

  f = open(InFileName, 'r')

  Z  = array('d')
  Bx = array('d')
  By = array('d')
  Bz = array('d')

  for l in f:
    [z, bx, by, bz] = map(float, l.split())
    z /= 1000.
    if z > ZMin and z < ZMax:
      Z.append(z)
      Bx.append(bx)
      By.append(by)
      Bz.append(bz)
      npZ += 1


  npXY = 1
  npZ  = len(Z)

  ZStep = (ZMax - ZMin) / npZ

    

  return SRWLMagFld3D(Bx, By, Bz, npXY, npXY, npZ, 0.0, 0.0, (npZ)*ZStep, 1, 1, None, None, _arZ=Z)





