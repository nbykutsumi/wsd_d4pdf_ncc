# %%
from numpy import *
from detect_fsub import *
from datetime import datetime, timedelta
import util
import numpy as np
import IO_Master
import ConstCyclone
import calendar
import detect_func
import collections
import os, sys
from importlib import import_module
detectName = (os.path.abspath(__file__)).split("/")[-2]
Cyclone = import_module("%s.Cyclone"%(detectName))

#-----------------------------------------
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

iDTime = datetime(2000,1,1,0)
#iDTime = datetime(2010,12,1,0)
eDTime = datetime(2010,12,31,18)

#iDTime = datetime(2006,1,1,6)   # HAPPI
#eDTime = datetime(2014,12,31,18)  # HAPPI

#print 'AAAAAAAAAAAAAAAAAAa'
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
  iDTime = datetime(iYear,iMon,1,0) # 2018/10/25
  eDTime = datetime(eYear,eMon,eDay,18)

#-------------------------

dDTime = timedelta(hours=6)
ret_lDTime = {False: util.ret_lDTime
             ,True : util.ret_lDTime_noleap
             }[noleap]

lDTime   = ret_lDTime(iDTime, eDTime, dDTime)

miss   = -9999.0
miss_int = -9999

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



#['lat','lon','nowpos','age','dura','epos','idate','ipos','nextpos','prepos','slp','slp_mean_adj','slp_mean_box','vortlw','vortlw_max_box','x','y','time']

ddtype={
    'lat':'float32',
    'lon':'float32',
    'nowpos':'int32',
    'age':'int32',
    'dura':'int32',
    'epos':'int32',
    'idate':'int32',
    'ipos':'int32',
    'nextpos':'int32',
    'prepos':'int32',
    'slp':'float32',
    'slp_mean_adj':'float32',
    'slp_mean_box':'float32',
    'vortlw':'float32',
    'vortlw_max_box':'float32',
    'x':'int32',
    'y':'int32',
    'time':'int32',
    'iedist':'float32',
}



#*********** functions *********
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass

def ret_a1iedist(a1ipos, a1epos):
  miss     = -9999.
  a1nx   = ones(len(a1ipos))*nx
  a1miss = ones(len(a1ipos))*miss_int
  lixy   = array(list(map(detect_func.fortpos2pyxy, a1ipos, a1nx, a1miss)))
  lexy   = array(list(map(detect_func.fortpos2pyxy, a1epos, a1nx, a1miss)))

  a1ix     = array(lixy[:,0], int32)
  a1iy     = array(lixy[:,1], int32)
  a1ex     = array(lexy[:,0], int32)
  a1ey     = array(lexy[:,1], int32)

  a1ilon   = a1lon.take(ma.masked_equal(a1ix, miss_int).filled(0))
  a1elon   = a1lon.take(ma.masked_equal(a1ex, miss_int).filled(0))
  a1ilat   = a1lat.take(ma.masked_equal(a1iy, miss_int).filled(0))
  a1elat   = a1lat.take(ma.masked_equal(a1ey, miss_int).filled(0))

  a1ilon   = ma.masked_where(a1ix==miss_int, a1ilon).filled(miss)
  a1elon   = ma.masked_where(a1ex==miss_int, a1elon).filled(miss)
  a1ilat   = ma.masked_where(a1iy==miss_int, a1ilat).filled(miss)
  a1elat   = ma.masked_where(a1ey==miss_int, a1elat).filled(miss)

  liedist  = list(map(detect_fsub.hubeny_real, a1ilat, a1ilon, a1elat, a1elon))# m
  a1iedist = array(liedist, float32)/1000.0
  a1iedist = ma.masked_invalid(a1iedist).filled(miss)
  return a1iedist

#********************************
a1lat    = iom.Lat
a1lon    = iom.Lon
a2lon, a2lat = meshgrid(a1lon, a1lat)
ny       = iom.ny
nx       = iom.nx
#---------
#--- a2pos ---------------------
# pos = 1,2,3,4,....
a2nowpos  = array(arange( ny*nx).reshape(ny,nx), int32) + 1

#-----------------------------------------
for idt, DTime in enumerate(lDTime):
  #******* init **********
  year = DTime.year
  mon  = DTime.month
  day  = DTime.day
  hour = DTime.hour
  if ((idt==0)or(DTime.month != lDTime[idt-1].month)):
    print(DTime.year, DTime.month)
    a2num  = zeros([ny,nx],float32).reshape(ny,nx)
    #--------------
    #lstype_ex  = ["dura","pgrad","vortlw","lat","lon","ipos","epos","idate","nowpos","time","prepos","nextpos"]
    lstype_ex  = ['lat','lon','nowpos','age','dura','epos','idate','ipos','nextpos','prepos','slp','slp_mean_adj','slp_mean_box','vortlw','vortlw_max_box','x','y','time']

    da1     = {}

    #for stype in lstype_ex + lstype_tc:
    for stype in lstype_ex:
      #if stype in ["dura","ipos","idate","nowpos","time","prepos","nextpos"]:
      #  da1[stype] = array([],int32  )
      #else:
      #  da1[stype] = array([],float32)
      da1[stype] = collections.deque([])
  #-------------
  #a2dura          = detect_fsub.solvelife_dura(a2life.T, miss_int).T 
  #---- shrink ---------------------
  a1age_tmp   = cy.load_clist_org('age',DTime) 
  a1dura_tmp  = cy.load_clist_org('dura',DTime) 
  a1idate_tmp = cy.load_clist_org('idate',DTime) 
  a1ipos_tmp         = cy.load_clist_org('ipos',DTime) 
  a1epos_tmp         = cy.load_clist_org('epos',DTime) 
  a1prepos_tmp       = cy.load_clist_org('prepos',DTime) 
  a1nextpos_tmp      = cy.load_clist_org('nextpos',DTime) 
  a1slp_tmp          = cy.load_clist_org('slp',DTime) 
  a1slp_mean_adj_tmp = cy.load_clist_org('slp_mean_adj',DTime) 
  a1slp_mean_box_tmp = cy.load_clist_org('slp_mean_box',DTime) 
  a1vortlw_tmp         = cy.load_clist_org('vortlw',DTime) 
  a1vortlw_max_box_tmp = cy.load_clist_org('vortlw_max_box',DTime) 
  a1x_tmp      = cy.load_clist_org('x',DTime) 
  a1y_tmp      = cy.load_clist_org('y',DTime) 
  a1lat_tmp    = a1lat[a1y_tmp]
  a1lon_tmp    = a1lon[a1x_tmp]
  a1nowpos_tmp = a2nowpos[a1y_tmp, a1x_tmp]

  #--- Screenig ----
  a1flag = ma.masked_greater_equal(a1dura_tmp, 0).mask

  a1age_tmp             = a1age_tmp           [a1flag]  
  a1dura_tmp            = a1dura_tmp          [a1flag]  
  a1idate_tmp           = a1idate_tmp         [a1flag] 
  a1ipos_tmp            = a1ipos_tmp          [a1flag]   
  a1epos_tmp            = a1epos_tmp          [a1flag]   
  a1prepos_tmp          = a1prepos_tmp        [a1flag]   
  a1nextpos_tmp         = a1nextpos_tmp       [a1flag]     
  a1slp_tmp             = a1slp_tmp           [a1flag]   
  a1slp_mean_adj_tmp    = a1slp_mean_adj_tmp  [a1flag]
  a1slp_mean_box_tmp    = a1slp_mean_box_tmp  [a1flag]
  a1vortlw_tmp          = a1vortlw_tmp        [a1flag] 
  a1vortlw_max_box_tmp  = a1vortlw_max_box_tmp[a1flag]
  a1x_tmp               = a1x_tmp             [a1flag]
  a1y_tmp               = a1y_tmp             [a1flag]
  a1lat_tmp             = a1lat_tmp           [a1flag]
  a1lon_tmp             = a1lon_tmp           [a1flag]
  a1nowpos_tmp          = a1nowpos_tmp        [a1flag]

  #--- a1time ------
  time            = year*10**6 + mon*10**4 + day*10**2 + hour
  a1time_tmp      = np.full(len(a1slp_tmp), time ,int32)

  #-----------------
  da1["age"    ].extend( a1age_tmp     )
  da1["dura"   ].extend( a1dura_tmp    )
  da1["idate"  ].extend( a1idate_tmp   )
  da1["ipos"   ].extend( a1ipos_tmp    )
  da1["epos"   ].extend( a1epos_tmp    )
  da1["prepos" ].extend( a1prepos_tmp    )
  da1["nextpos"].extend( a1nextpos_tmp    )
  da1["nowpos" ].extend( a1nowpos_tmp    )
  da1["slp"    ].extend( a1slp_tmp    )
  da1["slp_mean_adj"].extend( a1slp_mean_adj_tmp )
  da1["slp_mean_box"].extend( a1slp_mean_box_tmp )
  da1["vortlw"  ].extend( a1vortlw_tmp   )
  da1["vortlw_max_box"].extend( a1vortlw_max_box_tmp   )
  da1["x"     ].extend( a1x_tmp      )
  da1["y"     ].extend( a1y_tmp      )
  da1["lat"     ].extend( a1lat_tmp      )
  da1["lon"     ].extend( a1lon_tmp      )
  da1["time"    ].extend( a1time_tmp     )

  if ((DTime==lDTime[-1])or(DTime.month != lDTime[idt+1].month)):
    ##---- iedist ------------------------
    da1["iedist"] = ret_a1iedist(da1["ipos"], da1["epos"])
  
    #----- make dir ----
    sodir  = cy.path_clist("ipos", year, mon)[0]
    mk_dir(sodir)
  
    #---- save --
    _lstype = lstype_ex + ["iedist"]
    #_lstype = lstype_ex
    
    for stype in _lstype:
      datatype = ddtype[stype]
      soname = cy.path_clist(stype, year, mon)[1]
      a1out  = np.array(da1[stype]).astype(datatype)
      np.save(soname, a1out)
    print(soname)
  
