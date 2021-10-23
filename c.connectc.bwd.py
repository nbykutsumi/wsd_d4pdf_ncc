# %%
from numpy import *
from detect_fsub import *
from datetime import datetime, timedelta
import numpy as np
import util
import IO_Master
import ConstCyclone
import calendar
import os, sys
from importlib import import_module
detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))

##***************************
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

iDTime = datetime(2000,1,1,0)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2006,1,1,6)   # HAPPI
#eDTime = datetime(2014,12,31,18)   # HAPPI

#flgresume  = True
flgresume  = False
dbbaseDir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, noleap, flagresume, dbbaseDir, wsbaseDir = largv[1:1+8]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()

  if flagresume=="True": flagresume=True
  elif flagresume=="False": flagresume=False
  else: print("check flagresume",flagresume); sys.exit()

  iYear,iMon, eYear, eMon = list(map(int,largv[1+8:]))

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
lDTimeRev= lDTime[::-1]


#iom    = IO_Master.IO_Master(cfg, prj, model, run, res)
iom    = IO_Master.IO_Master(prj, model, run, res, dbbaseDir)

const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon
wsDir = wsbaseDir + '/%s'%(run)
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

#****************
miss_dbl     = -9999.0
miss_int     = -9999
endh         = 18
thdp         = 0.0  #[Pa]
thdist_search = 500.0*1000.0   #[m]
#####################################################
# functions
#####################################################
def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

#####################################################
def pyxy2fortpos(ix, iy, nx):
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
  if noleap==True:
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
#-- const --- 
thtopo   = cy.thtopo

#****************************************************
# read lat, lon data
#----------------------
a1lat, a1lon = iom.Lat, iom.Lon
ny           = iom.ny
nx           = iom.nx
X,Y          = meshgrid(arange(nx), arange(ny))
#**************************************************
# read topo data
a2topo      = iom.Load_const("topo")
a2mask_topo = ma.masked_greater(a2topo, thtopo)
#*************************************************
counter = 0
nlDTime = len(lDTimeRev)
for idt, DTime1 in enumerate(lDTimeRev):

  if ((DTime1==lDTimeRev[0])or(DTime1.month !=lDTime[idt-1])):
    print("connectc.bwd.py, backward",DTime1.year, DTime1.month)

  counter= counter + 1
  DTime0 = DTime1 - dDTime
  if ((noleap==True)&(DTime0.month==2)&(DTime0.day==29)):
    DTime0 = DTime0 - timedelta(days=1)

  #***************************************
  # "nextpos" for final timestep
  #**********
  if counter == 1:
    #------
    if flgresume == True:
      #preposnextname1 = cy.path_a2dat("prepos",DTime1+dDTime).srcPath
      #a2preposnext1   = fromfile(preposnextname1, int32).reshape(ny,nx)

      preposnextname1 = cy.path_clist_org("prepos",DTime1+dDTime)[1]
      ynextname1      = cy.path_clist_org("y",DTime1+dDTime)[1]
      xnextname1      = cy.path_clist_org("x",DTime1+dDTime)[1]

      a1preposnext1 = np.load(preposnextname1)
      a1ynext1      = np.load(ynextname1)
      a1xnext1      = np.load(xnextname1)

      a2preposnext1 = np.full([ny,nx],miss_int, int32)
      a2preposnext1[a1ynext1, a1xnext1] = a1preposnext1

      a2nextpos1      = np.full([ny,nx],miss_int, int32)
      for iynext in range(0, ny):
        for ixnext in range(0, nx):
          if (a2preposnext1[iynext, ixnext] != miss_int):
            (ix1,iy1) = fortpos2pyxy(a2preposnext1[iynext,ixnext], nx, miss_int)
            a2nextpos1[iy1,ix1] = pyxy2fortpos(ixnext, iynext, nx)
    #------
    else:
      a2nextpos1    = np.full([ny,nx], miss_int, int32)


    #nextposdir1   = cy.path_a2dat("nextpos",DTime1).srcDir
    #nextposname1  = cy.path_a2dat("nextpos",DTime1).srcPath

    yname1 = cy.path_clist_org("y",DTime1)[1]
    xname1 = cy.path_clist_org("x",DTime1)[1]
    a1y1   = np.load(yname1)
    a1x1   = np.load(xname1)

    nextposdir1  = cy.path_clist_org("nextpos",DTime1)[0]
    nextposname1 = cy.path_clist_org("nextpos",DTime1)[1]
    mk_dir(nextposdir1)
    #------
    a1nextpos1 = a2nextpos1[a1y1,a1x1]
    np.save(nextposname1, a1nextpos1)
    #a2nextpos1.tofile(nextposname1)
  #----------

  #***************************************
  #* names for 1
  #---------------------------------------
  nextposname1 = cy.path_clist_org("nextpos",DTime1)[1]
  preposname1  = cy.path_clist_org("prepos",DTime1)[1]
  agename1     = cy.path_clist_org("age",  DTime1)[1]
  yname1       = cy.path_clist_org("y",  DTime1)[1]
  xname1       = cy.path_clist_org("x",  DTime1)[1]

  #nextposname1 = cy.path_a2dat("nextpos",DTime1).srcPath
  #preposname1  = cy.path_a2dat("prepos",DTime1).srcPath
  #agename1     = cy.path_a2dat("age",  DTime1).srcPath
  #----------
  # read data
  #**********
  try:
    a1nextpos1   = np.load(nextposname1)
    a1prepos1    = np.load(preposname1)
    a1age1       = np.load(agename1)
    a1y1         = np.load(yname1)
    a1x1         = np.load(xname1)

    a2nextpos1   = np.full([ny,nx], miss_int, int32)
    a2prepos1    = np.full([ny,nx], miss_int, int32)
    a2age1       = np.full([ny,nx], miss_int, int32)

    a2nextpos1[a1y1,a1x1] = a1nextpos1
    a2prepos1 [a1y1,a1x1] = a1prepos1
    a2age1    [a1y1,a1x1] = a1age1

    #a2nextpos1   = fromfile(nextposname1,int32).reshape(ny,nx)
    #a2prepos1    = fromfile(preposname1, int32).reshape(ny,nx)
    #a2age1       = fromfile(agename1,    int32).reshape(ny,nx)
  except IOError:
    counter = counter -1
    print("No File:")
    print(preposname1)
    print(nextposname1)
    print(agename1)
    sys.exit()

  #**************************************
  #   inverse trace
  #--------------------------------------
  if (counter == 1):
    if flgresume == True:
      duraname2 = cy.path_clist_org("dura",DTime1+dDTime)[1]
      eposname2 = cy.path_clist_org("epos",DTime1+dDTime)[1]
      yname2    = cy.path_clist_org("y",DTime1+dDTime)[1]
      xname2    = cy.path_clist_org("x",DTime1+dDTime)[1]

      a1dura2   = np.load(duraname2)
      a1epos2   = np.load(eposname2)
      a1y2      = np.load(yname2)
      a1x2      = np.load(xname2)

      #duraname2 = cy.path_a2dat("dura",DTime1+dDTime).srcPath
      #a2dura2   = fromfile(duraname2, int32).reshape(ny,nx)

      #eposname2 = cy.path_a2dat("epos",DTime1+dDTime).srcPath
      #a2epos2   = fromfile(eposname2, int32).reshape(ny,nx)


      a2duranext= np.full([ny,nx],miss_int, int32)
      a2eposnext= np.full([ny,nx],miss_int, int32)
      for iy1 in range(ny):
        for ix1 in range(nx):
          if (a2nextpos1[iy1,ix1] !=miss_int): 
            ix2,iy2 = fortpos2pyxy(a2nextpos1[iy1,ix1], nx, miss_int)
            a2duranext[iy1,ix1] = a2dura2[iy2,ix2]
            a2eposnext[iy1,ix1] = a2epos2[iy2,ix2]
    else:  
      a2duranext   = np.full([ny,nx] ,miss_int, int32)
      a2eposnext   = np.full([ny,nx] ,miss_int, int32)

  #--------------------------
  # initialize a2dura1 and a2dura2_new
  #*****************
  a2dura1        = np.full([ny,nx], miss_int, int32) 
  a2duranext_new = np.full([ny,nx], miss_int, int32) 
  a2epos1        = np.full([ny,nx], miss_int, int32) 
  a2eposnext_new = np.full([ny,nx], miss_int, int32) 

  a2nextpos0     = np.full([ny,nx], miss_int, int32) 
  #*****************
  ax1 = ma.masked_where(a2age1 ==miss_int, X).compressed()
  ay1 = ma.masked_where(a2age1 ==miss_int, Y).compressed()
  for iy1,ix1 in zip(ay1, ax1):
    age1    = a2age1[iy1, ix1]
    duranext = a2duranext[iy1, ix1]
    eposnext = a2eposnext[iy1, ix1]
    (ix0,iy0) = fortpos2pyxy(a2prepos1[iy1,ix1], nx, miss_int)
    #---- 
    if (duranext == miss_int):
      #dura1 = 1000000* age1 + int(pgmax1)
      dura1 = age1
      epos1 = pyxy2fortpos(ix1,iy1,nx)
    else:
      dura1 = duranext
      epos1 = eposnext
    #----
    a2dura1[iy1, ix1] = dura1
    a2epos1[iy1, ix1] = epos1
    #-----------------------
    # fill a2duranext_new, a2eposnext_new
    #***************
    if (ix0 != miss_int):
      a2duranext_new[iy0, ix0] = dura1
      a2eposnext_new[iy0, ix0] = epos1
    #-----------------------
    # make "a2nextpos0"
    #***************
    if (iy0 != miss_int):
      a2nextpos0[iy0, ix0] = pyxy2fortpos(ix1, iy1, nx)

  #-------------------
  # replace a2duranext with new data
  #*******************
  a2duranext = a2duranext_new
  a2eposnext = a2eposnext_new
  #**************************************
  # write to file
  #--------------------------------------
  # out dir
  #**********
  duradir1     = cy.path_clist_org("dura"   ,DTime1)[0]
  eposdir1     = cy.path_clist_org("epos"   ,DTime1)[0]

  nextposdir0  = cy.path_clist_org("nextpos",DTime0)[0]

  mk_dir(duradir1)
  mk_dir(eposdir1)

  mk_dir(nextposdir0)
  #----------
  # out name
  #**********
  duraname1    = cy.path_clist_org("dura"   ,DTime1)[1]
  eposname1    = cy.path_clist_org("epos"   ,DTime1)[1]

  nextposname0 = cy.path_clist_org("nextpos",DTime0)[1]
  #----------
  # write to file
  #**********
  a1dura1 = a2dura1[a1y1, a1x1]
  a1epos1 = a2epos1[a1y1, a1x1]
  np.save(duraname1, a1dura1)
  np.save(eposname1, a1epos1)

  if idt != nlDTime-1:
    a1y0    = np.load(cy.path_clist_org("y" ,DTime0)[1])
    a1x0    = np.load(cy.path_clist_org("x" ,DTime0)[1])
    a1nextpos0 = a2nextpos0[a1y0, a1x0]
    np.save(nextposname0, a1nextpos0)

  #a2dura1.tofile(duraname1)
  #a2epos1.tofile(eposname1)
  #a2nextpos0.tofile(nextposname0) 

  print(duraname1)



# %%
