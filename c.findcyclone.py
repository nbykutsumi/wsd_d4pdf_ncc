from numpy import *
#from .detect_fsub import *
from detect_fsub import *
from datetime import datetime, timedelta
#from .detect_func import box_filtering
from detect_func import box_filtering
import util
import grids
import IO_Master
import ConstCyclone
#import Cyclone
import calendar
import os, sys, shutil
import numpy as np
import scipy.ndimage.filters as filters

import matplotlib
import matplotlib.pyplot as plt
from importlib import import_module

detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))
#--------------------------------------------------
#prj     = "JRA55"
#model   = "__"
#run     = "__"
#res     = "145x288"
#noleap  = False

#prj     = "HAPPI"
#model   = "MIROC5"
##run     = "C20-ALL-001-100"
#run     = "XXXX"
#res     = "128x256"
#noleap  = True

prj     = "d4PDF"
model   = "__"
#run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
run     = "XX-HPB-001"   # {expr}-{scen}-{ens}
res     = "320x640"
noleap  = False
dbbaseDir  = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
wsbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM'

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
  #iDTime = datetime(iYear,iMon,1,6)
  iDTime = datetime(iYear,iMon,1,0)
  eDTime = datetime(eYear,eMon,eDay,18)
#-------------------------


dDTime = timedelta(hours=6)

ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)

tstp        = "6hr"

#cfg          = config.cfg
#cfg['prj']   = prj    # for ConstCyclone
#cfg['model'] = model  # for ConstCyclone
#cfg['outbaseDir'] = cfg['baseDir'] + '/%s'%(run)
#iom    = IO_Master.IO_Master(cfg, prj, model, run, res)
wsDir = wsbaseDir + '/%s'%(run)
iom    = IO_Master.IO_Master(prj, model, run, res, dbbaseDir)

const  = ConstCyclone.Const(prj=prj, model=model)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon
cy     = Cyclone.Cyclone(baseDir=wsDir, const=const)

a1lat  = iom.Lat
a1lon  = iom.Lon
ny     = iom.ny
nx     = iom.nx
miss   = -9999.0

dlon = a1lon[1] - a1lon[0]
dlat = (a1lat[1:] - a1lat[:-1]).mean()
dx_4deg = int(4.0/dlon)
dy_4deg = int(4.0/dlat)
######################################################
#def box_filtering(a2in, a1y, a1x, dy, dx, func=None, cover_poles=True, miss_in=-9999., miss_out=-9999.):
#
#  if type(a2in) is np.ma.core.MaskedArray:
#    a2in = a2in.filled(miss_in)
#
#  a2large = grids.expand_map_global_2d(a2in, dy, dx, cover_poles=cover_poles)
#  a1y_large = a1y+dy
#  a1x_large = a1x+dx
#  ldyx = [[idy,idx] for idy in range(-dy,dy+1) for idx in range(-dx, dx+1)]
#  nybox = 2*dy+1
#  nxbox = 2*dx+1
#
#  a2container = np.empty([len(a1y), nybox*nxbox], a2in.dtype)
#  for i,(idy,idx) in enumerate(ldyx):
#    a2container[:, i] = a2large[a1y_large+idy, a1x_large+idx]
#
#  if func=='sum':
#    a1out = ma.masked_equal(a2container, miss_in).sum(axis=1)
#  elif func=='mean':
#    a1out = ma.masked_equal(a2container, miss_in).mean(axis=1)
#  elif func=='min':
#    a1out = ma.masked_equal(a2container, miss_in).min(axis=1)
#  elif func=='max':
#    a1out = ma.masked_equal(a2container, miss_in).max(axis=1)
#  else:
#    print 'check func',func
#    sys.exit()
#
#  return a1out.filled(miss_out)
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
#################################################
def mk_dir_tail(var, tstp, model, expr, ens):
  odir_tail = var + "/" + tstp + "/" +model + "/" + expr +"/"\
       +ens
  return odir_tail
#####################################################
def mk_namehead(var, tstp, model, expr, ens):
  namehead = var + "_" + tstp + "_" +model + "_" + expr +"_"\
       +ens
  return namehead
#****************************************************
def read_txtlist(iname):
  f = open(iname, "r")
  lines = f.readlines()
  f.close()
  lines = list(map(float, lines))
  aout  = array(lines, float32)
  return aout
#****************************************************
##**************************************************
## Mean Sea Level Pressure
##------------------------
##------------------------
for DTime in lDTime:
  year = DTime.year
  mon  = DTime.month
  day  = DTime.day
  hour = DTime.hour
  #pgraddir = cy.path_a2dat("pgrad",datetime(year,mon,1)).srcDir
  #mk_dir(pgraddir)
  #***************************************
  # pgrad
  #---------------------------------------
  #pgradname = cy.path_a2dat("pgrad",DTime).srcPath
  a2slp   = iom.Load_6hrSfc("slp", DTime)
  #findcyclone_out = detect_fsub.findcyclone_bn(a2slp.T, a1lat, a1lon, iom.miss,  miss)
  #a2pgrad = findcyclone_out.T
  #a2pgrad.tofile(pgradname)

  a2center = -detect_fsub.find_localmax(-a2slp.T, a1lat, a1lon, iom.miss,  0).T

  a1y,a1x = np.where( ma.masked_not_equal(a2center,0))

  a1slp = a2center[a1y,a1x]

  xdir, xname = cy.path_clist_org("x",DTime)
  ydir, yname = cy.path_clist_org("y",DTime)

  slpdir, slpname = cy.path_clist_org("slp",DTime)

  mk_dir(xdir)
  mk_dir(ydir)
  mk_dir(slpdir)
  np.save(xname, a1x)
  np.save(yname, a1y)
  np.save(slpname, a1slp)
  #print slpname
  #***************************************
  # Mean slp
  #---------------------------------------
  a1slp_mean_adj = box_filtering(a2slp, a1y, a1x, 1, 1, func='mean', miss_in=-9999, miss_out=-9999.)
  slpmean_adj_dir, slpmean_adj_name = cy.path_clist_org("slp_mean_adj",DTime)
  mk_dir(slpmean_adj_dir)
  np.save(slpmean_adj_name, a1slp_mean_adj)


  a1slp_mean_box = box_filtering(a2slp, a1y, a1x, dy_4deg, dx_4deg, func='mean', miss_in=-9999, miss_out=-9999.)
  slpmean_box_dir, slpmean_box_name = cy.path_clist_org("slp_mean_box",DTime)
  mk_dir(slpmean_box_dir)
  np.save(slpmean_box_name, a1slp_mean_box)

  ##***************************************
  ## rvort @ 850
  ##---------------------------------------
  rvortdir, rvortname = cy.path_clist_org("vortlw",DTime)
  mk_dir(rvortdir)

  a2u       = iom.Load_6hrPlev("ua", DTime, 850)
  a2v       = iom.Load_6hrPlev("va", DTime, 850)
  a2rvort   = detect_fsub.mk_a2rvort(a2u.T, a2v.T, a1lon, a1lat, iom.miss).T

  a2mask    = ma.masked_equal(a2rvort, iom.miss).mask
  a2rvort[:int(ny/2)] = -a2rvort[:int(ny/2)]  # The signs of the missing values in the south hemisphere are also fliped
  a2rvort   = ma.masked_where(a2mask, a2rvort).filled(iom.miss)

  #- find maximum vorticity in 9x9 box --
  a1maxvort_adj = box_filtering(a2rvort, a1y, a1x, 1, 1, func='max', miss_in=-9999, miss_out=-9999.)

  ##- find maximum vorticity in a (2*dy+1) x (2*dx+1) box --
  a1maxvort_box = box_filtering(a2rvort, a1y, a1x, dy_4deg, dx_4deg, func='max', miss_in=-9999, miss_out=-9999.)


  #-- Save ----
  vort_max_box_dir, vort_max_box_name = cy.path_clist_org("vortlw_max_box",DTime)
  vort_max_adj_dir, vort_max_adj_name = cy.path_clist_org("vortlw",DTime)


  mk_dir(vort_max_box_dir)
  mk_dir(vort_max_adj_dir)
  np.save(vort_max_box_name, a1maxvort_box)
  np.save(vort_max_adj_name, a1maxvort_adj)
  print(vort_max_box_name)

