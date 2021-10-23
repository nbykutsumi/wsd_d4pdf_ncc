# %%
from numpy import *
from   datetime import datetime, timedelta
import util
import config
import os, sys
import calendar
import IO_Master
import ConstCyclone
import Cyclone
import numpy as np
#******************************************************
#******************************************************
#prj     = "JRA55"
#model   = "__"
#run     = "__"
#res     = "145x288"
#tstp    = "6hr"
#noleap  = False

#prj     = "HAPPI"
#model   = "MIROC5"
##run     = "C20-ALL-001"
#run     = "XXXX"
#res     = "128x256"
#tstp    = "day"
#noleap  = True

prj     = "d4PDF"
model   = "__"
#run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
run     = "XX-HPB-001"   # {expr}-{scen}-{ens}
res     = "320x640"
tstp    = '6hr'
noleap  = False
dbbaseDir  = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
wsbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

iDTime = datetime(2010,1,1)
eDTime = datetime(2010,1,15)
#iDTime = datetime(2006,1,1,6)  # HAPPI
#eDTime = datetime(2016,1,1,0)  # HAPPI

#-- argv ----------------
largv = sys.argv
if len(largv)>1:
  prj, model, run, res, tstp, noleap, dbbaseDir, wsbaseDir = largv[1:1+8]
  if noleap=="True": noleap=True
  elif noleap=="False": noleap=False
  else: print "check noleap",noleap; sys.exit()

  iYear,iMon, eYear, eMon = map(int,largv[1+8:])

  eDay1 = calendar.monthrange(eYear,eMon)[1]
  iDTime = datetime(iYear,iMon,1,0)
  eDTime = datetime(eYear,eMon,eDay1,18)
#-------------------------

print "*"*50
print __file__[0]
print prj, model, run, res, tstp, noleap

lvar   = ["ua","va"]
plev   = 500

lhour  = {"6hr": [0,6,12,18]
         ,"day": [0]
         }[tstp]

dDTime = timedelta(days=1)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)

miss   = -9999.
#******************************************************
# set dlyrange
#******************************************************
dnx    = {}
dny    = {}
#****************************************************
#cfg    = config.cfg
#cfg['prj']   = prj    # for ConstCyclone
#cfg['model'] = model  # for ConstCyclone
#
#iom    = IO_Master.IO_Master(cfg, prj, model, run, res)
iom    = IO_Master.IO_Master(prj, model, run, res, dbbaseDir)

#const  = ConstCyclone.Const(cfg)
const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon
#cfg['outbaseDir'] = cfg['baseDir'] + '/%s'%(run)
wsDir = wsbaseDir + '/%s'%(run)

cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)
nx     = iom.nx
ny     = iom.ny

#dw         = 7
dw         = 3
ldaydelta  = range(-dw, dw+1)

if   tstp=="6hr": Load_Var = iom.Load_6hrPlev
elif tstp=="day": Load_Var = iom.Load_dayPlev
#####################################################
# Function
#####################################################
def check_file(sname):
  if not os.access(sname, os.F_OK):
    print "no file:",sname
    sys.exit()
#####################################################
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#******************************************************
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
        target = target + timedelta(days=-1)
      elif (target >= leapdate ) & (leapdate > today):
        target = target + timedelta(days=1)
  #-----------
  return target
  
#******************************************************
for var in lvar:
  #------
  #odir_root = os.path.join(cfg["outbaseDir"],"run.mean",var)
  odir_root = os.path.join(wsDir,"run.mean",var)
  #------------------------------
  # make heads and tails
  #------------------------------
  for DTime in lDTime:
    print DTime
    #*************
    year   = DTime.year
    mon    = DTime.month
    day    = DTime.day
    odir   = odir_root + "/%04d/%02d"%(year, mon)
    mk_dir(odir)
    #*************
    stime  = "%04d%02d%02d%02d"%(year,mon,day, 0)
    #***********
    oname  = odir + "/run.mean.%s.%04dhPa.%s.%s"%(var, plev, stime, res)
    #*********************
    # start running mean
    #*********************
    # container
    #********
    a3out  = zeros([len(ldaydelta)*len(lhour),ny,nx], float32)
    #********
    i=-1
    for daydelta in ldaydelta:
      target     = date_slide( year, mon, day, daydelta, noleap)
      targetyear = target.year
      targetmon  = target.month
      targetday  = target.day
      #-------------------
      for targethour in lhour:
        i = i+1
        tDTime = datetime(targetyear, targetmon, targetday, targethour)
        try:
          #ain    = iom.Load_6hrPlev(var, tDTime, plev)
          ain  = Load_Var(var, tDTime, plev)
        #except IOError:
        except:
          ain  = np.ones([ny,nx],'float32')*miss

        #--------------------
        # add 
        #--------------------
        a3out[i] = ain
    #*****************
    aout = ma.masked_equal(a3out, miss).mean(axis=0).astype('float32')
    if ma.isMA(aout):
      aout = aout.filled(miss)
    #*****************
    aout.tofile(oname)




# %%
