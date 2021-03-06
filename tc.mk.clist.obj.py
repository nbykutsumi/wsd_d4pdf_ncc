from numpy import *
from collections import deque
from datetime import datetime, timedelta
from detect_fsub import *
from detect_func import box_filtering
import detect_func
import calendar
import util
import IO_Master
import ConstCyclone
import numpy as np
import os, sys
from importlib import import_module
detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))

##--------------------------------------
#prj   = "JRA55"
#model = "__"
#run   = "__"
#res   = "145x288"
#plev_up  = 250
#noleap= False

#prj     = "HAPPI"
#model   = "MIROC5"
#run     = "C20-ALL-001"
#res     = "128x256"
#plev_up  = 250
#noleap  = True

prj     = "d4PDF"
model   = "__"
#run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
run     = "XX-HPB-001"   # {expr}-{scen}-{ens}
res     = "320x640"
plev_up  = 300
noleap  = False
dbbaseDir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
wsbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'


iDTime = datetime(2000,1,1,0)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2006,1,1,6)
#eDTime = datetime(2015,1,1,0)

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, noleap, dbbaseDir, wsbaseDir = largv[1:1+7]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()

  iYear,iMon, eYear, eMon = list(map(int,largv[1+7:]))
  eDay   = calendar.monthrange(eYear,eMon)[1]
  iDTime = datetime(iYear,iMon,1,0)   # 2020/04/28
  eDTime = datetime(eYear,eMon,eDay,18)
#-------------------------

dDTime = timedelta(hours=6)
ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)
lDTimeRev= lDTime[::-1]

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

#singleday = True
singleday = False
unitdist  = 10.0 # km / hour
#unitdist  = 150.0 # km / hour  # test
#----------------
miss  = -9999.
#lvar  = ["rvort", "dtlow", "dtmid", "dtup", "wmeanlow", "wmeanup", "sst","land"]
#lvar  = ["dtlow", "dtmid", "dtup", "wmeanlow", "wmeanup", "sst","land"]
lvar  = ["dtlw", "dtmd", "dtup", "wmaxlw", "wmaxup", "sst","land"]

a1lat = cy.Lat
a1lon = cy.Lon

dlon = a1lon[1] - a1lon[0]
dlat = (a1lat[1:] - a1lat[:-1]).mean()

dx_4deg = int(4.0/dlon)
dy_4deg = int(4.0/dlat)

def save_clist(var, a1dat, Year, Mon):
  clistDir, clistPath = cy.path_clist(var, Year, Mon)
  detect_func.mk_dir(clistDir)
  dNumType            = cy.dNumType

  array(a1dat, dtype=dNumType[var]).tofile(clistPath)  
  print(clistPath)

def ret_a2pgrad(DTime):
  return cy.load_a2dat("pgrad", DTime) 

#----------------------------
for idt, DTime in enumerate(lDTime):
  print(DTime)

  if ((idt==0)or(DTime.month != lDTime[idt-1].month)):
    Year = DTime.year
    Mon  = DTime.month
    #*** init ***********
    da1 = {}
    for var in lvar:
      da1[var]  = deque([])
    #********************
    # SST
    #-------------------- 
    a2sst = iom.Load_monSfc("sst", Year, Mon)  # t2m is called for d4PDF

    #********************
    # Land
    #-------------------- 
    a2land= iom.Load_const("land")
  #-------------------- 
  a1y    = cy.load_clist_org('y', DTime)
  a1x    = cy.load_clist_org('x', DTime)

  #-- SST and Land ----
  a1sst    = a2sst[a1y,a1x]
  a1land   = a2land[a1y,a1x]

  #-- Screening -------
  a1dura   = cy.load_clist_org('dura', DTime)
  a1flag   = ma.masked_greater_equal(a1dura,0).mask

  a1y    = a1y[a1flag]
  a1x    = a1x[a1flag]
  a1sst  = a1sst[a1flag]
  a1land = a1land[a1flag]

  #--------------------
  a2tlw   = iom.Load_6hrPlev("ta",DTime,850)
  a2tmd   = iom.Load_6hrPlev("ta",DTime,500)
  a2tup   = iom.Load_6hrPlev("ta",DTime,plev_up)
  a2ulw   = iom.Load_6hrPlev("ua",DTime,850)
  a2uup   = iom.Load_6hrPlev("ua",DTime,plev_up)
  a2vlw   = iom.Load_6hrPlev("va",DTime,850)
  a2vup   = iom.Load_6hrPlev("va",DTime,plev_up)


  
  #-- Maximum and mean temperature ----------
  a1tmaxup = box_filtering(a2tup, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss, miss_out=miss) 

  a1tmaxmd = box_filtering(a2tmd, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss, miss_out=miss) 

  a1tmaxlw = box_filtering(a2tlw, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss, miss_out=miss) 

  a1taveup = box_filtering(a2tup, a1y, a1x, dy_4deg, dx_4deg, func='mean', miss_in=miss, miss_out=miss) 

  a1tavemd = box_filtering(a2tmd, a1y, a1x, dy_4deg, dx_4deg, func='mean', miss_in=miss, miss_out=miss) 

  a1tavelw = box_filtering(a2tlw, a1y, a1x, dy_4deg, dx_4deg, func='mean', miss_in=miss, miss_out=miss) 

  a1dtup   = (ma.masked_equal(a1tmaxup,miss) - ma.masked_equal(a1taveup,miss)).filled(miss)

  a1dtmd   = (ma.masked_equal(a1tmaxmd,miss) - ma.masked_equal(a1tavemd,miss)).filled(miss)

  a1dtlw   = (ma.masked_equal(a1tmaxlw,miss) - ma.masked_equal(a1tavelw,miss)).filled(miss)

  #-- Maximum wind ----------
  a2wup   = np.sqrt( np.square(ma.masked_equal(a2uup,miss)) + np.square(ma.masked_equal(a2vup,miss)) )

  a2wlw   = np.sqrt( np.square(ma.masked_equal(a2ulw,miss)) + np.square(ma.masked_equal(a2vlw,miss)) )

  a1wmaxup = box_filtering(a2wup, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss, miss_out=miss) 

  a1wmaxlw = box_filtering(a2wlw, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=miss, miss_out=miss) 

  #--------------------
  da1["dtlw"  ].extend( a1dtlw  )
  da1["dtmd"  ].extend( a1dtmd  )
  da1["dtup"  ].extend( a1dtup  )
  da1["wmaxlw"].extend( a1wmaxlw)
  da1["wmaxup"].extend( a1wmaxup)
  da1["sst"   ].extend( a1sst   )
  da1["land"  ].extend( a1land  )


  #- write clist --
  if ((DTime==lDTime[-1])or(DTime.month != lDTime[idt+1].month)):
    for var in lvar:
      a1out = np.array(da1[var], float32)
      outpath = cy.path_clist(var, Year, Mon)[1]
      np.save(outpath, a1out)
      print(len(a1out), a1out.dtype, outpath)
