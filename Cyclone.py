from numpy import *
import numpy as np
#from .detect_fsub import *
from detect_fsub import *
from datetime import datetime, timedelta
from collections import deque, OrderedDict
import os
import itertools
import socket
#import BestTrackTC
from . import BestTrackTC
import calendar
#****************************************************
def ret_latlonbnd(a1lat, a1lon):
  dlat0 = a1lat[1] -a1lat[0]
  dlat1 = a1lat[-1]-a1lat[-2]
  dlon0 = a1lon[1] -a1lon[0]
  dlon1 = a1lon[-1]-a1lon[-2]

  ilat = np.array([a1lat[0]])
  elat = np.array([a1lat[-1]])
  ilon = np.array([a1lon[0]])
  elon = np.array([a1lon[-1]])

  a1latbnd = 0.5*(a1lat[:-1] + a1lat[1:])
  a1lonbnd = 0.5*(a1lon[:-1] + a1lon[1:])

  a1latbnd = np.concatenate([ilat-dlat0*0.5, a1latbnd, elat+dlat1*0.5])
  a1lonbnd = np.concatenate([ilon-dlon0*0.5, a1lonbnd, elon+dlon1*0.5])

  return a1latbnd, a1lonbnd


def read_txtlist(iname):
  f = open(iname, "r")
  lines = f.readlines()
  f.close()
  lines = list(map(float, lines))
  aout  = array(lines, float32)
  return aout

def ret_lYM(iYM, eYM):
  iYear, iMon = iYM
  eYear, eMon = eYM
  if iYear==eYear:
    lYM = [[iYear,Mon] for Mon in range(iMon,eMon+1)]
  else:
    lYM = [[iYear,Mon] for Mon in range(iMon,12+1)]
  
    lYM = lYM + [[Year,Mon] \
                  for Year in range(iYear+1,eYear)\
                  for Mon  in range(1,12+1)]

    lYM = lYM + [[eYear,Mon] for Mon in range(1,eMon+1)]
  return lYM

def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

def ret_lDTime_noleap(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )

  TstpLeap = []
  for Year in range(iDTime.year, eDTime.year+1):
      if calendar.isleap(Year):
          itstp = ((datetime(Year,2,29,0)-iDTime).total_seconds()-1)\
                 /dDTime.total_seconds() + 1
          etstp = ((datetime(Year,3,1,0)-iDTime).total_seconds()-1)\
                 /dDTime.total_seconds()

          if itstp < 0: itstp=0
          if etstp > total_steps-1: etstp = total_steps-1

          TstpLeap = TstpLeap + list(range(int(itstp), int(etstp)+1))

  ltstp = [i for i in range(total_steps) if i not in TstpLeap]

  return [iDTime + dDTime*i for i in ltstp]

def ret_lDTime_fromYM(iYM,eYM, dDTime, hour0, noleap=False):
  iYear,iMon = iYM
  eYear,eMon = eYM

  eDay = calendar.monthrange(eYear,eMon)[1]
  iDTime = datetime(iYear,iMon,1,hour0)
  eDTime = datetime(eYear,eMon,eDay,23)
  if noleap is False:
      return ret_lDTime(iDTime, eDTime, dDTime)
  else:
      return ret_lDTime_noleap(iDTime, eDTime, dDTime)


#---------------------------------------------------
def solve_time(stime):
  year = int( stime/10**6 )
  mon  = int( (stime - year*10**6)/10**4 )
  day  = int( (stime - year*10**6 - mon*10**4)/10**2)
  hour = int( (stime - year*10**6 - mon*10**4 - day*10**2) )
  return year, mon, day, hour
#---------------------------------------------------
def fortpos2pyxy(number, nx, miss_int):
  if (number == miss_int):
    iy_py = miss_int
    ix_py = miss_int
  else:
    iy_py = int((number-1.0)/nx)      # iy_py = 0,1,2,..
    ix_py = int(number - nx*iy_py-1)   # ix_py = 0,1,2,..
  #----
  return ix_py, iy_py
#---------------------------------------------------
def expand8grids(lxy, ny, nx):
  # lxy = [(x1,y1), (x2,y2), ...]
  if len(lxy) == 0:
    lXY = array([])
  else:
    lXY = deque([])
    for x,y in lxy:
      lXY.extend(list(itertools.product([x-1,x,x+1],[y-1,y,y+1])))
  
    ax, ay = list(zip(*lXY))
    ax = ma.masked_equal(ax,-1).filled(nx-1)
    ax = ma.masked_equal(ax,nx).filled(0)
    ay = ma.masked_equal(ay,-1).filled(0)
    ay = ma.masked_equal(ay,ny).filled(ny-1)
  
    lXY = array(list(zip(ax,ay)))
  return  lXY
  
#---------------------------------------------------
class Cyclone():
  def __init__(self, baseDir, const):
    self.baseDir = baseDir
    basebaseDir  = "/".join(self.baseDir.split("/")[:-1])

    self.thtopo = const['thtopo']
    self.thdura = const['thdura']
    self.thsst  = const['thsst']
    #self.thpgrad= const['thpgrad']
    self.exrvort= const['exrvort']
    self.tcrvort= const['tcrvort']
    self.thwcore= const['thwcore']
    self.thwind = const['thwind']
    self.thwdif = const['thwdif']
    #------------

    self.Lat     = const['Lat'] 
    self.Lon     = const['Lon'] 

    self.ny      = len(self.Lat)
    self.nx      = len(self.Lon)
    self.dNumType= {"life"    :int32,
                    "dura"    :int32,
                    "ipos"    :int32,
                    "epos"    :int32,
                    "idate"   :int32,
                    "time"    :int32,
                    "age"     :int32,
                    "nowpos"  :int32,
                    "lastpos" :int32,
                    "prepos"  :int32,
                    "nextpos" :int32,
                    "pgrad"   :float32,
                    "vortlw"  :float32,
                    "pgmax"   :float32,
                    "pmean"   :float32,
                    #"iedist"  :float32,
                    #"rvort"   :float32,
                    "dtlow"   :float32,
                    "dtmid"   :float32,
                    "dtup"    :float32,
                    "wmeanlow":float32,
                    "wmeanup" :float32,
                    "sst"     :float32,
                    "land"    :float32,
                    "initsst" :float32,
                    "initland":float32,
                    "iedist"  :float32
                   }



  def path_clist_org(self, var, DTime):
    Year = DTime.year
    Mon  = DTime.month
    Day  = DTime.day
    Hour = DTime.hour
    tstep   = "6hr"
    srcDir  = os.path.join(self.baseDir, tstep, var, "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04d%02d%02d%02d.%dx%d.npy"%(var,Year,Mon,Day,Hour, self.ny, self.nx))
    return srcDir, srcPath

  def load_clist_org(self, var, DTime):
    srcPath = self.path_clist_org(var, DTime)[1]
    return np.load(srcPath)

  def load_a2dat(self, var, DTime):
    srcPath = self.path_a2dat(var, DTime)
    a2dat = fromfile(self.srcPath, self.dNumType[var]).reshape(self.ny,self.nx)
    return a2dat


  def path_a2dat(self, var, DTime):
    Year = DTime.year
    Mon  = DTime.month
    Day  = DTime.day
    Hour = DTime.hour
    self.tstep   = "6hr"
    self.srcDir  = os.path.join(self.baseDir, self.tstep, var, "%04d"%(Year), "%02d"%(Mon))
    self.srcPath = os.path.join(self.srcDir, "%s.%04d%02d%02d%02d.%dx%d"%(var,Year,Mon,Day,Hour, self.ny, self.nx))
    return self

  def load_a2dat(self, var, DTime):
    srcPath = self.path_a2dat(var, DTime)
    a2dat = fromfile(self.srcPath, self.dNumType[var]).reshape(self.ny,self.nx)
    return a2dat

  def path_clist(self, var, Year, Mon):
    self.tstep   = "6hr"
    srcDir  = os.path.join(self.baseDir, self.tstep, "clist", "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04d.%02d.npy"%(var,Year,Mon))
    return srcDir, srcPath

  def load_clist(self, var, Year, Mon):
    srcPath = self.path_clist(var, Year, Mon)[1]
    return np.load(srcPath)

  #def dictC(self, Year, Mon, varname="pgrad", tctype="obj"):
  # if tctype == "obj":
  #   return self.dictC_objTC(Year, Mon, varname=varname)
  # if tctype == "bst":
  #   return self.dictC_bstTC(Year, Mon, varname=varname)

  def mkInstDictC(self, iYM, eYM, varname="vortlw", tctype="bst"):
   if tctype == "obj":
     return self.mkInstDictC_objTC(iYM, eYM, varname=varname)
   if tctype == "bst":
     return self.mkInstDictC_bstTC(iYM, eYM, varname=varname)
   if tctype == "notc":
     return self.mkInstDictC_noTC(iYM, eYM, varname=varname)



  #def dictC_bstTC(self, Year, Mon, varname="pgrad"):
  #  thpgrad   = self.thpgrad
  #  thdura    = self.thdura

  #  dictExC   = {}
  #  dictTC    = {}
  #  da1       = {}

  #  #--- TC dictionary ----
  #  bst       = BestTrackTC.BestTrack("IBTrACS")
  #  dictTC    = bst.ret_dpyxy(Year, self.Lon, self.Lat)

 
  #  #---- make ExC dictionary -----
  #  lvar      = ["dura","pgrad","nowpos","time"]
  #  for var in lvar:
  #     da1[var]  = self.load_clist(var, Year, Mon)

  #  nlist    = len(da1["dura"])

  #  for i in range(nlist):
  #    dura        = da1["dura"    ][i]
  #    pgrad       = da1["pgrad"   ][i]
  #    nowpos      = da1["nowpos"  ][i]
  #    time        = da1["time"    ][i]
  #    oVar        = da1[varname][i]
  #    #---- dura -------
  #    if dura < thdura:
  #      #print "dura",dura,"<",thdura
  #      continue
  #    #---- thpgrad ----
  #    if pgrad < thpgrad:
  #      #print "pgrad",pgrad,"<",thpgrad
  #      continue
  #
  #    #---- time -------
  #    Year,Mon,Day,Hour = solve_time(time)
  #    DTime             = datetime(Year,Mon,Day,Hour)
  #
  #    #---- nowpos  ----
  #    x,y               = fortpos2pyxy(nowpos, self.nx, -9999)


  #    #-- check TC ----
  #    lxyTC  = expand8grids(dictTC[DTime], self.ny, self.nx)
  #    if [x,y] in lxyTC:
  #      continue 
 
  #    #-- List  -------
  #    oList  = [x,y,oVar] 


  #    try:
  #      dictExC[DTime].append(oList)
  #    except KeyError:
  #      dictExC[DTime] = [oList]
  #  #--------------------------------------------
  #  self.dictTC  = dictTC
  #  self.dictExC = dictExC
  #  return self


  def mkInstDictC_noTC(self, iYM, eYM, varname="vortlw"):
    thdura    = self.thdura

    dictExC   = OrderedDict()
    dictTC    = OrderedDict()
    da1       = OrderedDict()

    iYear,iMon= iYM
    eYear,eMon= eYM
    lYM       = ret_lYM(iYM, eYM)

    #---- make ExC dictionary -----
    lvar      = ["dura","nowpos","time","vortlw"]
    if varname not in lvar: lvar = lvar + [varname] 

    for Year,Mon in lYM:
      for var in lvar:
         da1[var]  = self.load_clist(var, Year, Mon)
  
      nlist    = len(da1["dura"])
  
      for i in range(nlist):
        dura        = da1["dura"    ][i]
        nowpos      = da1["nowpos"  ][i]
        time        = da1["time"    ][i]
        oVar        = da1[varname][i]
        #---- dura -------
        if dura < thdura:
          #print "dura",dura,"<",thdura
          continue
   
        #---- time -------
        Year,Mon,Day,Hour = solve_time(time)
        DTime             = datetime(Year,Mon,Day,Hour)
    
        #---- nowpos  ----
        x,y               = fortpos2pyxy(nowpos, self.nx, -9999)
  
        #-- List  -------
        oList  = [x,y,oVar] 
  
        try:
          dictExC[DTime].append(oList)
        except KeyError:
          dictExC[DTime] = [oList]
    #--------------------------------------------
    return dictExC, dictTC



  def mkInstDictC_bstTC(self, iYM, eYM, varname="vortlw"):
    thdura    = self.thdura

    dictExC   = OrderedDict()
    dictTC    = OrderedDict()
    da1       = OrderedDict()

    iYear,iMon= iYM
    eYear,eMon= eYM
    lYM       = ret_lYM(iYM, eYM)

    #a1lat = self.Lat
    #a1lon = self.Lon
    #a1latbnd, a1lonbnd = ret_latlonbnd(a1lat, a1lon)

    #--- TC dictionary ----
    bst       = BestTrackTC.BestTrack("IBTrACS")
    for Year in range(iYear,eYear+1):
      iDTime = datetime(Year,1,1,0)
      eDTime = datetime(Year,12,31,18)


      dictTC.update(bst.ret_dpyxy(iDTime, eDTime, self.Lon, self.Lat))

    #---- make ExC dictionary -----
    lvar      = ["dura","nowpos","time","vortlw"]
    if varname not in lvar: lvar = lvar + [varname] 

    for Year,Mon in lYM:
      for var in lvar:
         da1[var]  = self.load_clist(var, Year, Mon)
  
      nlist    = len(da1["dura"])
  
      for i in range(nlist):
        dura        = da1["dura"    ][i]
        nowpos      = da1["nowpos"  ][i]
        time        = da1["time"    ][i]
        oVar        = da1[varname][i]
        #---- dura -------
        if dura < thdura:
          #print "dura",dura,"<",thdura
          continue
   
        #---- time -------
        Year,Mon,Day,Hour = solve_time(time)
        DTime             = datetime(Year,Mon,Day,Hour)
    
        #---- nowpos  ----
        x,y               = fortpos2pyxy(nowpos, self.nx, -9999)
  
        #-- check TC ----
        lxyTC  = expand8grids(dictTC[DTime], self.ny, self.nx)
        if [x,y] in lxyTC:
          continue 
   
        #-- List  -------
        oList  = [x,y,oVar] 
  
        try:
          dictExC[DTime].append(oList)
        except KeyError:
          dictExC[DTime] = [oList]
    #--------------------------------------------
    return dictExC, dictTC


  def mkInstDictC_objTC(self, iYM, eYM, varname="vortlw"):
    thrvort   = self.tcrvort
    thwcore   = self.thwcore
    thdura    = self.thdura
    thinitsst = self.thsst 
    thwind    = self.thwind
    thwdif    = self.thwdif  # wind_low - wind_up

    dictExC   = OrderedDict()
    dictTC    = OrderedDict()
    da1       = OrderedDict()

    #---- Initialize ----
    ''' No-leap is not considered '''

    lDTime = ret_lDTime_fromYM(iYM,eYM,timedelta(hours=6), hour0=0, noleap=False)

    for DTime in lDTime:
      dictExC[DTime] = []
      dictTC[DTime]  = []

    #--------------------

    #lvar      = ["dura","pgrad","nowpos","nextpos","time","iedist","rvort","dtlow","dtmid","dtup","initsst","initland"]
    #lvar      = ["dura","nowpos","time","vortlw","dtlow","dtmid","dtup","initsst","initland"]
    lvar      = ["x","y","dura","time","vortlw","dtlw","dtmd","dtup","initsst","initland","wmaxup","wmaxlw"] 
   
    if varname not in lvar: lvar = lvar + [varname] 

    lYM       = ret_lYM(iYM, eYM)
    for Year,Mon in lYM:
      for var in lvar:
        if var=='dtsum':
          da1[var]  = self.load_clist('dtlw', Year, Mon)\
                     +self.load_clist('dtmd', Year, Mon)\
                     +self.load_clist('dtup' , Year, Mon)

        else:
          da1[var]  = self.load_clist(var, Year, Mon)

      #---- make dictionary -----
      nlist    = len(da1["dura"])
      for i in range(nlist):
        y           = da1['y'][i]
        x           = da1['x'][i]
        dura        = da1["dura"    ][i]
        #nowpos      = da1["nowpos"  ][i]
        time        = da1["time"    ][i]
        #iedist      = da1["iedist"  ][i]
        #rvort       = abs(da1["vortlw"   ][i])
        rvort       = da1["vortlw"   ][i]
        dtlow       = da1["dtlw"   ][i]
        dtmid       = da1["dtmd"   ][i]
        dtup        = da1["dtup"    ][i]
        initsst     = da1["initsst" ][i]
        initland    = da1["initland"][i]
        #nextpos     = da1["nextpos" ][i]
        wup         = da1["wmaxup" ][i]
        wlow        = da1["wmaxlw" ][i] 

        oVar        = da1[varname][i]
        #---- dura -------
        if dura < thdura:
          #print "dura",dura,"<",thdura
          continue
    
        #---- time -------
        Year,Mon,Day,Hour = solve_time(time)
        DTime             = datetime(Year,Mon,Day,Hour)
    
        #---- nowpos  ----
        #x,y               = fortpos2pyxy(nowpos, self.nx, -9999)
    
        #-- List  -------
        oList  = [x,y,oVar] 
        #---- thrvort ----
        if rvort < thrvort:
          #print "rvort",rvort,"<",thrvort
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]
    
          continue
    
        #---- thwcore ----
        if dtup + dtmid + dtlow < thwcore:
        #if (dtup <thwcore)or(dtmid<thwcore)or(dtmid<thwcore):
          #print "thwcore",dtup+dtmid+dtlow,"<",thwcore
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]
    
          continue

        #---- maximum wind @850 --
        if wlow < thwind:
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]

          #print "wup > wlow !!"
          continue 

        #---- wup & wlow --
        if (wlow - wup) < thwdif:
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]

          #print "wup > wlow !!"
          continue 
 
        #print 'low,up',wlow,wup, wlow-wup
        #---- initsst ----
        if initsst < thinitsst:
          #print "initsst",initsst,"<",thinitsst
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]
    
          continue
    
        #---- initland ----
        if initland > 0.0:
          #print "initland",initland,">",0.0
          try:
            dictExC[DTime].append(oList)
          except KeyError:
            dictExC[DTime] = [oList]
    
          continue
    
        #---- TC ----
        try:
          dictTC[DTime].append(oList)
        except KeyError:
          dictTC[DTime] = [oList]
  
  
    #--------------------------------------------
    return dictExC, dictTC


class Cyclone_2D(Cyclone):
  def __init__(self, iYM, eYM, cfg, const, tctype, miss=-9999.):
    """
    cfg : Dictionary with "prj", "model", "run", "res"
    tctype = "bst" or "obj"
    """
    Cyclone.__init__(self, cfg, const)
    self.instDict  = self.mkInstDictC(iYM, eYM, varname="vortlw",tctype=tctype)
    self.a2miss    = ones([self.ny, self.nx], float32)*miss
    self.miss      = miss
    self.tctype    = tctype
    self.res       = cfg["res"]

  def mk_a2tc(self, DTime, locfill=False):
    if len(self.instDict.dictTC[DTime]) >0:
      aList     = list(zip(*(self.instDict.dictTC[DTime])))
      a2dat     = self.a2miss.copy()
      if locfill==False:
        a2dat[ aList[1], aList[0]] = aList[-1]
      else:
        a2dat[ aList[1], aList[0]] = locfill

    else:
      a2dat     = self.a2miss.copy()
    return a2dat

  def mk_a2exc(self, DTime, locfill=False):
    if len(self.instDict.dictExC[DTime]) >0:
      aList     = list(zip(*(self.instDict.dictExC[DTime])))
      a2dat     = self.a2miss.copy()

      if locfill==False:
        a2dat[ aList[1], aList[0]] = aList[-1]
      else:
        a2dat[ aList[1], aList[0]] = locfill

    else:
      a2dat     = a2miss.copy()
    return a2dat

  def mkMask_exc(self, DTime, radkm=1000, miss=False):
    if type(miss) == bool: miss=self.miss

    a2loc     = self.mk_a2exc(DTime)
    return detect_fsub.mk_territory(a2loc.T, self.Lon, self.Lat, radkm*1000., imiss=self.miss, omiss=miss).T

  def mkMask_tc(self, DTime, radkm=1000, miss=False):
    if type(miss) == bool: miss=self.miss
    a2loc     = self.mk_a2tc(DTime)
    return detect_fsub.mk_territory(a2loc.T, self.Lon, self.Lat, radkm*1000., imiss=self.miss, omiss=miss).T

  def path_exc(self, DTime):
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour
    srcDir  = os.path.join(self.baseDir, self.tstep, "exc.%stc"%(self.tctype), "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04d%02d%02d%02d.%s"%("exc",Year,Mon,Day,Hour,self.res))

    return srcDir, srcPath

  def path_tc(self, DTime):
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour

    srcDir  = os.path.join(self.baseDir, self.tstep, "tc.%s"%(self.tctype), "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04d%02d%02d%02d.%s"%("tc",Year,Mon,Day,Hour,self.res))
    return srcDir, srcPath

  def path_Mask_exc(self, DTime, radkm=1000):
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour

    srcDir  = os.path.join(self.baseDir, self.tstep, "mask.exc.%stc"%(self.tctype), "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04dkm.%04d%02d%02d%02d.%s"%("exc",radkm,Year,Mon,Day,Hour,self.res))
    return srcDir, srcPath

  def path_Mask_tc(self, DTime, radkm=1000):
    Year    = DTime.year
    Mon     = DTime.month
    Day     = DTime.day
    Hour    = DTime.hour

    srcDir  = os.path.join(self.baseDir, self.tstep, "mask.tc.%s"%(self.tctype), "%04d"%(Year), "%02d"%(Mon))
    srcPath = os.path.join(srcDir, "%s.%04dkm.%04d%02d%02d%02d.%s"%("tc",radkm,Year,Mon,Day,Hour,self.res))
    return srcDir, srcPath


