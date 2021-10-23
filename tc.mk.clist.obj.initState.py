from numpy import *
from collections import deque
from   datetime import datetime, timedelta
import sys, os
import util
import IO_Master
import ConstCyclone
import numpy as np
from importlib import import_module
detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))


#prj   = "JRA55"
#model = "__"
#run   = "__"
#res   = "125x288"
#noleap= False
##iYear_data = 1958
##iMon_data  = 1
#iYear_data = 2004
#iMon_data  = 1

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

iYear_data = 2000
iMon_data  = 1
iYear, iMon = [2000,1]
eYear, eMon = [2010,12]

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, noleap, dbbaseDir, wsbaseDir = largv[1:1+7]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print("check noleap",noleap); sys.exit()

  iYear,iMon, eYear, eMon = list(map(int,largv[1+7:1+7+4]))
  iYear_data, iMon_data   = list(map(int,largv[1+7+4:1+7+4+2]))
#-------------------------
iom    = IO_Master.IO_Master(prj, model, run, res, dbbaseDir)

const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon

wsDir = wsbaseDir + '/%s'%(run)
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)
lYM = util.ret_lYM([iYear,iMon],[eYear,eMon])

#**********************************************
def ret_a1initState(var, year,mon,dinitState_pre):
  #----------
  dinitState              = {}
  dinitState[-9999,-9999] = -9999.0
  #----------
  a1idate     = cy.load_clist("idate",year,mon) 
  a1ipos      = cy.load_clist("ipos" ,year,mon) 
  a1time      = cy.load_clist("time" ,year,mon) 
  a1state     = cy.load_clist(var    ,year,mon) 
  a1land      = cy.load_clist("land" ,year,mon) 

  #------------------------
  print(len(a1idate), len(a1ipos), len(a1time), len(a1state), len(a1land))

  n  = len(a1idate)
  ldat    = deque([])
  for i in range(n):
    idate = a1idate[i]
    ipos  = a1ipos [i]
    time  = a1time [i]
    state = a1state[i]
    #print idate, ipos, time
    #--- check initial time --
    if time == idate:
      dinitState[idate, ipos] = state

    #-----------
    try:
      ldat.append( dinitState[idate, ipos] )
    except KeyError:
      try:
        ldat.append( dinitState_pre[idate, ipos])
      except:
        ldat.append( -9999.0)
#        sys.exit() 
  #---------------------------
  a1initState = array(ldat, float32)
  return dinitState, a1initState
#**********************************************

#--- init ----
date_first = datetime(iYear,iMon, 1)
date_pre   = date_first + timedelta(days = -2)
year_pre   = date_pre.year
mon_pre    = date_pre.month
if (iYear == iYear_data)&(iMon ==iMon_data):
  dinitsst   = {} 
  dinitland  = {} 
else:
  dinitsst , a1temp = ret_a1initState("sst" , year_pre, mon_pre, {} )
  dinitland, a1temp = ret_a1initState("land", year_pre, mon_pre, {} )
#-------------
for [year, mon] in lYM:
  dinitsst_pre          = dinitsst
  dinitsst, a1initsst   = ret_a1initState( "sst", year, mon, dinitsst_pre )

  dinitland_pre         = dinitland
  dinitland, a1initland = ret_a1initState( "land",year, mon, dinitland_pre )
 
 
  #---- oname ----------------
  name_sst  = cy.path_clist("initsst" ,year,mon)[1]
  name_land = cy.path_clist("initland",year,mon)[1]
  np.save(name_sst,  a1initsst)
  np.save(name_land, a1initland)
  #a1initsst.tofile(name_sst)
  #a1initland.tofile(name_land)
  print(name_sst)

  
