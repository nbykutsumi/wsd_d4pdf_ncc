from numpy import *
import numpy as np
from datetime import datetime, timedelta
import os, sys
import calendar
import socket
from bisect import bisect_right
#**********************************************
#def nearest_idx(aSrc,val):
#    ''' return nearest index. by HJKIM'''
#    if hasattr(val,'__iter__'): return [abs(aSrc-v).argmin() for v in val]
#    else: return abs(aSrc-val).argmin()

def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

#**********************************************
class IBTrACS(object):
  def __init__(self):
    #-- check host --
    hostname = socket.gethostname()
    if hostname == "well":
      #self.baseDir  = "/media/disk2/data/ibtracs"
      self.baseDir  = "/home/utsumi/mnt/lab_tank/utsumi/data/ibtracs"
    if hostname in ["shui","mizu","naam"]:
      self.baseDir  = "/tank/utsumi/data/ibtracs"
    #----------------
    #self.Versions= ["v03r04","v03r06","v03r08","v03r10"]
    self.Versions= ["v04r00"]


  def ret_path(self, ver="v04r00"):
    srcDir  = os.path.join(self.baseDir, ver)
    #self.srcPath = self.srcDir + "/Year.%04d.ibtracs_all.%s.csv"%(Year,ver)
    return srcDir + "/ibtracs.since1980.list.%s.csv"%(ver)

  def ret_dlonlat(self,iDTime,eDTime, ver="v04r00",dhours=6):
    #-- open Y=Year-1 ----
    srcPath = self.ret_path(ver=ver)
    if os.path.exists(srcPath):
      f       = open(srcPath, "r")
      lines   = f.readlines()[2:]
      f.close()
    else:
      print 'No file'
      print srcPath
      sys.exit()
    #--- init dict ---
    lDTime = ret_lDTime(iDTime,eDTime,timedelta(hours=dhours))
    dout   = {}
    for DTimeTmp in lDTime:
      dout[DTimeTmp] = []

    #-----------------
    for line in lines:
      line     = line.split(",")
      isotime  = line[6].split(" ")
      date     = map(int, isotime[0].split("-"))
      Year     = date[0]
      Mon      = date[1]
      Day      = date[2]
      Hour     = int(isotime[1].split(":")[0])
      DTimeTmp = datetime(Year,Mon,Day,Hour)
      #--- check Time --
      if (DTimeTmp<iDTime) or (eDTime<DTimeTmp):
        continue
      #--- check Hour --
      if Hour % dhours !=0: continue

      #--- check nature --
      nature   = line[7].strip()
      if nature not in ["TS"]:
        continue

      #-----------------
      tcname   = line[5].strip()
      tcid     = line[0]
      lat      = float(line[8])
      lon      = float(line[9])
      if (lon < 0.0):
        lon = 360.0 + lon
      #-----------------
      DTime    = datetime(Year,Mon,Day,Hour)
      dout[DTime].append([lon,lat])
    #---
    return dout
  #################################################
  def ret_dpyxy(self, iDTime, eDTime, a1lonbnd, a1latbnd, ver="v04r00"):
    ''' a1lonbnd, a1latbnd: boundaries of the grids, increasing order '''
    dlonlat = self.ret_dlonlat(iDTime,eDTime,ver)
    lkey    = dlonlat.keys()
    dout    = {}
    nx,ny   = len(a1lonbnd)-1, len(a1latbnd)-1

    for key in lkey:
      llonlat = dlonlat[key]
      if len(llonlat)==0:
        dout[key] = []
      else:
        aLon, aLat= zip(*llonlat)
        a1x = np.array([bisect_right(a1lonbnd, lon) for lon in aLon]) -1
        a1y = np.array([bisect_right(a1latbnd, lat) for lat in aLat]) -1
        a1x = ma.masked_outside(a1x, 0, nx-1).filled(-9999)
        a1y = ma.masked_outside(a1y, 0, ny-1).filled(-9999)
        dout[key] = zip( a1x, a1y )

    #---
    return dout


class IBTrACS_2D(IBTrACS):
  def __init__(self, iDTime,eDTime, a1lonbnd, a1latbnd, miss=-9999.0, ver="v04r00"):
    IBTrACS.__init__(self)
    self.dpyxy   = self.ret_dpyxy(iDTime, eDTime, a1lonbnd, a1latbnd, ver="v04r00")
    self.LatBnd  = a1latbnd # Grid boundaries
    self.LonBnd  = a1lonbnd # Grid boundaries, 0 - 360
    self.Lat     = (a1latbnd[:-1] + a1latbnd[1:])*0.5   # Grid centers 
    self.Lon     = (a1lonbnd[:-1] + a1lonbnd[1:])*0.5   # Grid centers
    self.ny      = len(self.Lat)
    self.nx      = len(self.Lon)
    self.a2miss  = ones([self.ny, self.nx], float32)*(miss)
    self.a2zero  = zeros([self.ny, self.nx], int32)

  def load_a2dat(self, DTime):
    lpyxy        = self.dpyxy[DTime]
    a2dat        = self.a2zero.copy()

    a1x, a1y = zip(*lpyxy)
    a1flag = np.logical_and(a1x>=0, a1y>=0)
    a1x = np.array(a1x)[a1flag]
    a1y = np.array(a1y)[a1flag] 
    a2dat[ a1y, a1x ] = 1
    return a2dat

    



