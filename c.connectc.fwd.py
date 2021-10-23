# %%
from numpy import *
import numpy as np
from detect_fsub import *
from datetime import datetime, timedelta
import util
#import config
import IO_Master
import ConstCyclone
import calendar
import os, sys
from importlib import import_module
detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))

##***************************
##--------------------------------------------------
#prj     = "JRA55"
#model   = "__"
#run     = "__"
#res     = "145x288"
#noleap  = False

#prj     = "HAPPI"
#model   = "MIROC5"
#run     = "C20-ALL-001"
#res     = "128x256"
#noleap  = True

prj     = "d4PDF"
model   = "__"
#run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
run     = "XX-HPB-001"   # {expr}-{scen}-{ens}
res     = "320x640"
noleap  = False
dbbaseDir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

#iDTime = datetime(2000,1,1,0)
iDTime = datetime(2010,1,31,0)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2006,1,1,6)  # HAPPI
#eDTime = datetime(2014,12,31,18)  # HAPPI

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, noleap, dbbaseDir, wsbaseDir = largv[1:1+7]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()
  iYear,iMon, eYear, eMon = list(map(int,largv[1+7:]))

  eDay   = calendar.monthrange(eYear,eMon)[1]
  #iDTime = datetime(iYear,iMon,1,6)
  iDTime = datetime(iYear,iMon,1,0)  # 2018/10/25
  eDTime = datetime(eYear,eMon,eDay,18)
#-------------------------

dDTime = timedelta(hours=6)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)


#cfg          = config.cfg
#cfg['prj']   = prj    # for ConstCyclone
#cfg['model'] = model  # for ConstCyclone
#cfg['outbaseDir'] = cfg['baseDir'] + '/%s'%(run)
#iom    = IO_Master.IO_Master(cfg, prj, model, run, res)
iom    = IO_Master.IO_Master(prj, model, run, res, dbbaseDir)

const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon
wsDir = wsbaseDir + '/%s'%(run)
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

#----------------
hinc         = 6
miss_dbl     = -9999.0
miss_int     = -9999
endh         = 18
#thpgrad      = const['thpgrad'] #[Pa/1000km]
#exrvort      = const['exrvort'] #[s-1]
#thpgrad      = const['thpgrad_min'] #[Pa/1000km]  2020/5/5
thpdif       = const['thpdif_min'] # Pa, 8deg box average - center
exrvort      = const['rvort_min']  #[s-1]  2020/5/5
thdist_search = 500.0*1000.0   #[m]
thtopo       = const['thtopo']

#####################################################
# functions
#####################################################
#####################################################
def fortxy2fortpos(ix, iy, nx):
  ix     = ix + 1  # ix = 1,2,.. nx
  iy     = iy + 1  # iy = 1,2,.. ny
  #number = iy* nx + ix +1
  number = (iy-1)* nx + ix
  return number
#####################################################
def fortpos2pyxy(number, nx, miss_int):
  if (number == miss_int):
    iy0 = miss_int
    ix0 = miss_int
  else:
    iy0 = int((number-1)/nx)  +1  # iy0 = 1, 2, ..
    ix0 = number - nx*(iy0-1)     # ix0 = 1, 2, ..

    iy0 = iy0 -1    # iy0 = 0, 1, .. ny-1
    ix0 = ix0 -1    # ix0 = 0, 1, .. nx-1
  #----
  return ix0, iy0
#####################################################
def check_file(sname):
  if not os.access(sname, os.F_OK):
    print("no file:",sname)
    sys.exit()
#####################################################
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#####################################################
def date_slide(year,mon,day, daydelta, noleap):
  today       = datetime(year, mon, day)
  target      = today + timedelta(daydelta)
  targetyear  = target.year
  #***********
  if noleap == True:
    if ( calendar.isleap(targetyear) ):
      leapdate   = datetime(targetyear, 2, 29)
      #---------
      if (target <= leapdate) & (leapdate < today):
        target = target + timedelta(-1)
      elif (target >= leapdate ) & (leapdate > today):
        target = target + timedelta(1)
  #-----------
  return target
#****************************************************
def read_txtlist(iname):
  f = open(iname, "r")
  lines = f.readlines()
  f.close()
  lines = list(map(float, lines))
  aout  = array(lines, float32)
  return aout
#******************************************************
#****************************************************
# read lat, lon data
#----------------------
a1lat, a1lon = iom.Lat, iom.Lon
ny           = iom.ny
nx           = iom.nx
#**************************************************
# read topo data
a2topo      = iom.Load_const("topo")
a2mask_topo = ma.masked_greater(a2topo, thtopo)
#**************************************************
counter = 0
for idt, DTime in enumerate(lDTime):
  if ((idt==0)or(DTime.month != lDTime[idt-1].month)):
    print("connectc.fwd.py", model, "forward",DTime.year, DTime.month)

  #***********************
  counter = counter + 1
  #---------
  DTime1 = DTime
  year1, mon1, day1, hour1 = DTime1.year, DTime1.month, DTime1.day, DTime1.hour

  DTime0 = DTime1 - timedelta(hours=hinc)

  if ((noleap==True)&(DTime0.month==2)&(DTime0.day==29)):
    DTime0 = DTime0 - timedelta(days=1)

  year0, mon0, day0, hour0 = DTime0.year, DTime0.month, DTime0.day, DTime0.hour
  stimed0= "%04d%02d%02d%02d"%(year0,mon0,day0,0)
 

  print(DTime, DTime0) 
  #***************************************
  #* names for 0
  #---------------------------------------
  slpname0     = cy.path_clist_org("slp",  DTime0)[1]
  pboxname0    = cy.path_clist_org("slp_mean_box",  DTime0)[1]
  locname0     = cy.path_clist_org("vortlw",  DTime0)[1]
  xname0       = cy.path_clist_org("x",  DTime0)[1]
  yname0       = cy.path_clist_org("y",  DTime0)[1]

  preposname0  = cy.path_clist_org("prepos",DTime0)[1]
  iposname0    = cy.path_clist_org("ipos",   DTime0)[1]
  idatename0   = cy.path_clist_org("idate",  DTime0)[1]
  agename0     = cy.path_clist_org("age",   DTime0)[1]
  #wdifname0    = cy.path_clist_org("wdif",   DTime0)[1]  # temp
  #timename0    = cy.path_clist_org("time",   DTime0)[1]  # temp

  uadir0       = os.path.join(cy.baseDir,"run.mean","ua","%02d"%(year0),"%02d"%(mon0))
  vadir0       = os.path.join(cy.baseDir,"run.mean","va","%02d"%(year0),"%02d"%(mon0))

  uaname0     =  os.path.join(uadir0, "run.mean.%s.0500hPa.%s.%s"%("ua", stimed0, res))
  vaname0     =  os.path.join(vadir0, "run.mean.%s.0500hPa.%s.%s"%("va", stimed0, res))

  #***************************************
  #* names for 1
  #---------------------------------------
  preposdir1  = cy.path_clist_org("prepos", DTime1)[0]
  iposdir1    = cy.path_clist_org("ipos",   DTime1)[0]
  idatedir1   = cy.path_clist_org("idate",  DTime1)[0]
  agedir1     = cy.path_clist_org("age",    DTime1)[0]
  #wdifdir1    = cy.path_clist_org("wdif",   DTime1)[0]

  mk_dir(preposdir1)
  mk_dir(iposdir1)
  mk_dir(idatedir1)
  mk_dir(agedir1)
  #mk_dir(wdifdir1)

  slpname1     = cy.path_clist_org("slp",    DTime1)[1]
  pboxname1    = cy.path_clist_org("slp_mean_box", DTime1)[1]
  locname1     = cy.path_clist_org("vortlw", DTime1)[1]
  xname1       = cy.path_clist_org("x", DTime1)[1]
  yname1       = cy.path_clist_org("y", DTime1)[1]

  preposname1  = cy.path_clist_org("prepos", DTime1)[1]
  iposname1    = cy.path_clist_org("ipos",   DTime1)[1]
  idatename1   = cy.path_clist_org("idate",  DTime1)[1]
  agename1     = cy.path_clist_org("age",   DTime1)[1]
  #wdifname1    = cy.path_clist_org("wdif",   DTime1)[1]

  #***************************************
  # read data
  #---------------------------------------
  #   for 0
  #************
  print("iposname0")
  print(iposname0)
  if ( os.access(iposname0, os.F_OK) ):
    a1slp0   = np.load(slpname0)
    a1pbox0  = np.load(pboxname0)
    a1loc0   = np.load(locname0)
    a1y0     = np.load(yname0)
    a1x0     = np.load(xname0)

    a1pdif0  = a1pbox0 - a1slp0

    a2pdif0  = np.full([ny,nx], miss_dbl, 'float32')
    a2loc0  = np.full([ny,nx], miss_dbl, 'float32')

    a2pdif0[a1y0,a1x0] = a1pdif0
    a2loc0 [a1y0,a1x0] = a1loc0

    a2ua0      = fromfile(uaname0,      float32).reshape(ny, nx)
    a2va0      = fromfile(vaname0,      float32).reshape(ny, nx)

    #--
    print(preposname0)
    a1prepos0 = np.load(preposname0)
    a1ipos0   = np.load(iposname0)
    a1idate0  = np.load(idatename0) 

    a2prepos0 = np.full([ny,nx], miss_int, 'int32')
    a2ipos0   = np.full([ny,nx], miss_int, 'int32')
    a2idate0  = np.full([ny,nx], miss_int, 'int32')

    a2prepos0[a1y0,a1x0] = a1prepos0
    a2ipos0  [a1y0,a1x0] = a1ipos0
    a2idate0 [a1y0,a1x0] = a1idate0

    #-- 
    #try: 
    #  a1age0    = np.load(agename0)
    #except IOError:
    #  a1age0    = np.load(timename0)
 
    a1age0    = np.load(agename0)
    a2age0    = np.full([ny,nx], miss_int, 'int32')
    a2age0[a1y0,a1x0] = a1age0

    #a1wdif0   = np.load(wdifname0)
    #a2wdif0   = np.full([ny,nx], miss_dbl, 'float32')
    #a2wdif0[a1y0,a1x0] = a1wdif0


  elif ( counter == 1):
    a2pdif0    = np.full([ny,nx], miss_dbl, float32)
    a2loc0     = np.full([ny,nx], miss_dbl, float32)
    a2ua0      = np.full([ny,nx], miss_dbl, float32)
    a2va0      = np.full([ny,nx], miss_dbl, float32)
    #--
    a2prepos0  = np.full([ny,nx], miss_int, int32)
    a2ipos0    = np.full([ny,nx], miss_int, int32)
    a2idate0   = np.full([ny,nx], miss_int, int32)
    a2age0     = np.full([ny,nx], miss_int, int32)

    #a2wdif0    = np.full([ny,nx], miss_dbl, float32)

  else:
    print("nofiles: DTime0 = ",DTime0)
    print("iposname0 =", iposname0)
    sys.exit()

  #if ( os.access(iposname0, os.F_OK) ):
  #  a2pgrad0   = fromfile(pgradname0,   float32).reshape(ny, nx)
  #  a2loc0     = fromfile(locname0,     float32).reshape(ny, nx)
  #  a2ua0      = fromfile(uaname0,      float32).reshape(ny, nx)
  #  a2va0      = fromfile(vaname0,      float32).reshape(ny, nx)
  #  #--
  #  a2prepos0 = fromfile(preposname0, int32).reshape(ny, nx)
  #  a2ipos0    = fromfile(iposname0,    int32).reshape(ny, nx)
  #  a2idate0   = fromfile(idatename0,   int32).reshape(ny, nx)
 
  #  #-- 
  #  try: 
  #    a2age0    = fromfile(agename0,    int32).reshape(ny, nx)
  #  except IOError:
  #    a2age0    = fromfile(timename0,   int32).reshape(ny, nx)
 
  #elif ( counter == 1):
  #  a2pgrad0    = array(ones(ny*nx).reshape(ny,nx)*miss_dbl, float32)
  #  a2loc0     = array(ones(ny*nx).reshape(ny,nx)*miss_dbl, float32)
  #  a2ua0      = array(zeros(ny*nx).reshape(ny,nx), float32)
  #  a2va0      = array(zeros(ny*nx).reshape(ny,nx), float32)
  #  #--
  #  a2prepos0 = array(ones(ny*nx).reshape(ny,nx)*miss_int, int32)
  #  a2ipos0    = array(ones(ny*nx).reshape(ny,nx)*miss_int, int32)
  #  a2idate0   = array(ones(ny*nx).reshape(ny,nx)*miss_int, int32)
  #  a2age0    = array(ones(ny*nx).reshape(ny,nx)*miss_int, int32)
  #else:
  #  print "nofiles: DTime0 = ",DTime0
  #  print "iposname0 =", iposname0
  #  sys.exit()
  #------------
  #   for 1
  #************
  a1slp1    = np.load(slpname1)
  a1pbox1   = np.load(pboxname1)
  a1loc1    = np.load(locname1)
  a1y1      = np.load(yname1)
  a1x1      = np.load(xname1)

  a1pdif1   = a1pbox1 - a1slp1

  a2pdif1   = np.full([ny,nx], miss_dbl, float32)
  a2loc1    = np.full([ny,nx], miss_dbl, float32)
  a2pdif1[a1y1,a1x1] = a1pdif1
  a2loc1 [a1y1,a1x1] = a1loc1

  #a2pgrad1     = fromfile(pgradname1, float32).reshape(ny, nx)
  #a2loc1       = fromfile(locname1, float32).reshape(ny, nx)
  #*********************
  # mask pdif < thpdif and high altitudes
  #*********************
  a2loc0   = ma.masked_where(a2pdif0 < thpdif, a2loc0)
  a2loc1   = ma.masked_where(a2pdif1 < thpdif, a2loc1)

  a2loc0   = ma.masked_where(a2mask_topo.mask, a2loc0).filled(miss_dbl)
  a2loc1   = ma.masked_where(a2mask_topo.mask, a2loc1).filled(miss_dbl)


  #****************************************
  # connectc
  ##***************************************
  #ctrackout = detect_fsub.connectc_bn_nopgmax(\
  #   a2pgrad0.T, a2pgrad1.T, a2ua0.T, a2va0.T\
  #   , a2ipos0.T, a2idate0.T, a2age0.T\
  #   , a1lon, a1lat, thdp, thdist_search, hinc, miss_dbl, miss_int\
  #   , year1, mon1, day1, hour1\
  #   )

  #ctrackout = detect_fsub.connectc_vort_inertia(\
  #   a2pgrad0.T, a2pgrad1.T, a2ua0.T, a2va0.T\
  #   , a2prepos0.T, a2ipos0.T, a2idate0.T, a2age0.T\
  #   , a1lon, a1lat, exrvort, thdist_search, hinc, miss_dbl, miss_int\
  #   , year1, mon1, day1, hour1\
  #   )

  ctrackout = detect_fsub.connectc_vort_inertia(\
     a2loc0.T, a2loc1.T, a2ua0.T, a2va0.T\
     , a2prepos0.T, a2ipos0.T, a2idate0.T, a2age0.T\
     , a1lon, a1lat, exrvort, thdist_search, hinc, miss_dbl, miss_int\
     , year1, mon1, day1, hour1\
     )

  a2prepos1 = array(ctrackout[0].T, int32)
  a2ipos1    = array(ctrackout[1].T, int32)
  a2idate1   = array(ctrackout[2].T, int32)
  a2age1    = array(ctrackout[3].T, int32)

  #****************************************
  # write to file
  #----------------------------------------
  a1prepos1 = a2prepos1[a1y1,a1x1]
  a1ipos1   = a2ipos1  [a1y1,a1x1]
  a1idate1  = a2idate1 [a1y1,a1x1]
  a1age1    = a2age1   [a1y1,a1x1]

  np.save(preposname1,  a1prepos1)
  np.save(iposname1,    a1ipos1)
  np.save(idatename1,   a1idate1)
  np.save(agename1,     a1age1)

  #a2prepos1.tofile(preposname1)
  #a2ipos1.tofile(iposname1)
  #a2idate1.tofile(idatename1)
  #a2age1.tofile(agename1)
  



# %%
