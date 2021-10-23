from datetime import datetime, timedelta
from calendar import timegm
from numpy    import ma, array, iterable, shape
import os, math
import calendar
import numpy as np

def dtime2tstmp(DTime):
  if iterable(DTime) ==True:
    Shape  = shape(DTime)
    return array([timegm(dtime.timetuple())
                  for dtime in DTime]).reshape(Shape)
  else:
    return timegm(DTime.timetuple())

def tstmp2dtime(Tstmp):
  if iterable(Tstmp) ==True:
    Shape  = shape(Tstmp)
    return array([datetime.utcfromtimestamp(tstmp)
                  for tstmp in Tstmp]).reshape(Shape)
  else:
    return datetime.utcfromtimestamp(Tstamp)

def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass

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

          TstpLeap = TstpLeap + range(int(itstp), int(etstp)+1)

  ltstp = [i for i in range(total_steps) if i not in TstpLeap]

  return [iDTime + dDTime*i for i in ltstp]


def ret_day_of_year(DTime):
    Year   = DTime.year
    DTime0 = datetime(Year,1,1,0)
    dDTime = DTime - DTime0
    return dDTime.days + 1


def ret_dMonName():
  return {1:"Jan",2:"Feb",3:"Mar",4 :"Apr",5 :"May",6 :"Jun"
         ,7:"Jul",8:"Aug",9:"Sep",10:"Oct",11:"Nov",12:"Dec"}

def ret_lYM(iYM, eYM):
  """
  iYM = [iYear, iMon], eYM = [eYear, eMon]
  """
  iYear, iMon = iYM
  eYear, eMon = eYM
  lYM = []
  for Year in range(iYear, eYear+1):
    if iYear == eYear:
      lMon = range(iMon,eMon+1)
    elif Year == iYear:
      lMon = range(iMon,12+1)
    elif Year == eYear:
      lMon = range(1,eMon+1)
    else:
      lMon = range(1,12+1)   

    for Mon in lMon:
      lYM.append([Year,Mon])
  return lYM 

def shift_YM(Year,Mon, dMon):
  # dMon: can be positive and negative

  Mon2 = Mon + dMon
  if (Mon2)%12 ==0:
    oMon   = 12
    oYear  = int(Year + math.floor(Mon2/12)-1)
  else:
    oMon   = (Mon + dMon)%12
    oYear  = int(Year + math.floor(Mon2/12))
  return oYear,oMon 


def array2csv(a):
  if len(a.shape)==1:
    sout = "\n".join(map(str,a))
  elif len(a.shape)==2:
    lline = [",".join( map(str,line)) for line in a]
    sout  = "\n".join(lline).strip()
  return sout

def array2tab(a):
  if len(a.shape)==1:
    sout = "\n".join(map(str,a))
  elif len(a.shape)==2:
    lline = ["\t".join( map(str,line)) for line in a]
    sout  = "\n".join(lline).strip()
  return sout

def list2csv(a):
  if type(a[0]) !=list:
    sout = "\n".join(map(str,a))
  elif type(a[0]) ==list:
    lline = [",".join( map(str,line)) for line in a]
    sout  = "\n".join(lline).strip()
  return sout

def list2tab(a):
  if type(a[0]) !=list:
    sout = "\n".join(map(str,a))
  elif type(a[0]) ==list:
    lline = ["\t".join( map(str,line)) for line in a]
    sout  = "\n".join(lline).strip()
  return sout


def join_list_cols(ll):
    lout = []
    for l in ll:
        lout.append(list(l))
    return map(list, zip(*lout))
         

def ret_lmon(season):
  if type(season)==int:
    lmon  = [season]
  elif type(season)==np.int64:
    lmon  = [season]
  elif season.isdigit():
    lmon  = [int(season)]
  elif season == "DJF":
    lmon  = [1,2, 12]
  elif season == "MAM":
    lmon  = [3,4,5]
  elif season == "JJA":
    lmon  = [6,7,8]
  elif season == "SON":
    lmon  = [9,10,11]
  elif season == "ALL":
    lmon  = [1,2,3,4,5,6,7,8,9,10,11,12]
  elif type(season) == int:
    lmon  = [season]
  elif season == "NDJFMA":
    lmon  = [11,12,1,2,3,4]
  elif season == "MJJASO":
    lmon  = [5,6,7,8,9,10]
  elif season == "JJASON":
    lmon  = [6,7,8,9,10,11]
  elif season == "JJAS":
    lmon  = [6,7,8,9]
  elif season == "JASO":
    lmon  = [7,8,9,10]
  elif season == "JFMA":
    lmon  = [1,2,3,4]
  elif season == "DJFM":
    lmon  = [12,1,2,3]
  elif season == "NoJune":
    lmon  = [1,2,3,4,5,7,8,9,10,11,12]
  elif season == "NoJJ":
    lmon  = [1,2,3,4,5,8,9,10,11,12]
  else:
    print "check season",season, type(season)
  return lmon


def nhourly2monthly(a3in, nh, calc="mean",miss_in=None, miss_out=None):
    # nh : n-hourly
    # calc: "mean","sum"
    # if miss_in !=None and miss_out==None, output will be filled with miss_in
    #
    dtype   = a3in.dtype
    per_day = int(24/nh)
    nstep, ny, nx = a3in.shape
    if nstep/per_day == 365:
        Year = 1999  # just for calc
    elif nstep/per_day == 366:
        Year = 2000  # ust for calc

    if miss_in !=None:
        a3in = ma.masked_equal(a3in, miss_in)

    a3mon = np.empty([12,ny,nx], dtype=dtype)
    for Mon in range(1,12+1):
        iDay = 1
        eDay = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,iDay,0)
        eDTime = datetime(Year,Mon,eDay,0)
        iDOY   = ret_day_of_year(iDTime)
        eDOY   = ret_day_of_year(eDTime)
        ik     = (iDOY-1)*per_day
        ek     = (eDOY-1+1)*per_day
        if calc=="mean":
            a2mon = a3in[ik:ek,:,:].mean(axis=0)
        elif calc=="sum":
            a2mon = a3in[ik:ek,:,:].sum(axis=0)
        elif calc=="std":
            a2mon = a3in[ik:ek,:,:].std(axis=0)
        else:
            print 'check calc', calc

        if ma.isMaskedArray(a2mon):
            if miss_out !=None:
                a2mon = a2mon.filled(miss_out)

        a3mon[Mon-1,:,:] = a2mon
    return a3mon


def nhourly2day2monthly(a3in, nh, calc="mean",miss_in=None, miss_out=None):
    # nh : n-hourly
    # calc: "mean","sum"
    # if miss_in !=None and miss_out==None, output will be filled with miss_in
    #
    dtype   = a3in.dtype
    per_day = int(24/nh)
    nstep, ny, nx = a3in.shape
    if nstep/per_day == 365:
        year = 1999  # just for calc
    elif nstep/per_day == 366:
        year = 2000  # just for calc

    if miss_in !=None:
        a3in = ma.masked_equal(a3in, miss_in)

    a3mon = np.empty([12,ny,nx], dtype=dtype)
    for mon in range(1,12+1):
        iday = 1
        eday = calendar.monthrange(year,mon)[1]
        idtime = datetime(year,mon,iday,0)
        edtime = datetime(year,mon,eday,0)
        idoy   = ret_day_of_year(idtime)
        edoy   = ret_day_of_year(edtime)
        ik     = (idoy-1)
        ek     = (edoy-1+1)
        a4tmp  = a3in.reshape(-1,per_day,ny,nx)



        if calc=="mean":
            a2mon = a4tmp[ik:ek,:,:].mean(axis=1).mean(axis=0)
        elif calc=="sum":
            a2mon = a4tmp[ik:ek,:,:].mean(axis=1).sum(axis=0)
        elif calc=="std":
            a2mon = a4tmp[ik:ek,:,:].mean(axis=1).std(axis=0)
        else:
            print 'check calc', calc

        if ma.isMaskedArray(a2mon):
            if miss_out !=None:
                a2mon = a2mon.filled(miss_out)

        a3mon[mon-1,:,:] = a2mon
    return a3mon

