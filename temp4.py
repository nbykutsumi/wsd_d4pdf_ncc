import matplotlib
matplotlib.use('Agg')
#----------------------------------
import sys, os
myPath  = os.path.abspath(__file__)
myDir   = os.path.dirname(myPath)
parDir  = "/".join(myDir.split("/")[:-1])
sys.path.append(parDir)
#----------------------------------
from   mpl_toolkits.basemap import Basemap
from   numpy import *
from   detect_fsub import *
from   datetime import datetime, timedelta
import matplotlib.pyplot as plt
import config
import detect_func
import IO_Master
import BestTrackTC
import util
import Cyclone
import ConstCyclone
#--------------------------------------
#prj     = "HAPPI"
#model   = "MIROC5"
#run     = "C20-ALL-001"
#res     = "128x256"

prj     = "d4PDF"
model   = "__"
run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
res     = "320x640"
noleap  = False

iDTime = datetime(2010,1,1,0)
eDTime = datetime(2010,1,5,0)

cfg          = config.cfg
cfg['prj']   = prj    # for ConstCyclone
cfg['model'] = model  # for ConstCyclone
cfg['outbaseDir'] = cfg['baseDir'] + '/%s'%(run)
iom    = IO_Master.IO_Master(cfg, prj, model, run, res)

const  = ConstCyclone.Const(cfg)
const['Lat'] = iom.Lat
const['Lon'] = iom.Lon
cy     = Cyclone.Cyclone(cfg, const)


#----------------------------------
a1lat   = iom.Lat
a1lon   = iom.Lon
ny      = len(a1lat)
nx      = len(a1lon)

lonlatfontsize = 10.0
#lonrotation    = 90
lonrotation    = 0
miss_int= -9999

[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
#------------------------
iYM = [iDTime.year,iDTime.month]
eYM = [eDTime.year,eDTime.month]
_, dtcxy  = cy.mkInstDictC_objTC(iYM,eYM,varname='vortlw')
_, dtcppos= cy.mkInstDictC_objTC(iYM,eYM,varname='prepos')
#------------------------
print "Basemap"
figmap   = plt.figure()
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

for DTime in lDTime:
  lxyz = dtcxy[DTime]
  lprepos = dtcppos[DTime]

  if len(lxyz)==0: continue

  for i,(x,y,z) in enumerate(lxyz):
    lon = a1lon[x]
    lat = a1lat[y]

    xpre,ypre = Cyclone.fortpos2pyxy(lprepos[i][2], nx, miss_int=-9999)
    if (xpre>0):
      lonpre,latpre = a1lon[xpre], a1lat[ypre]
    else:
      lonpre,latpre = lon, lat

    #if ((lllon<lon)&(lon<urlon)&(lllat<lat)&(lat<urlat)):
    #  print x,y, "***",lon, lat
    #  M.plot( lon, lat, "o")

    #------------------------------------
    lon1, lat1 = lon, lat
    lon2, lat2 = lonpre, latpre

    scol = 'r'
    if abs(lon1 - lon2) >= 180.0:
      #--------------
      if (lon1 > lon2):
        lon05_1  = 360.0
        lon05_2  = 0.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      elif (lon1 < lon2):
        lon05_1  = 0.0
        lon05_2  = 360.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      #--------------
      M.plot( (lon1, lon05_1), (lat1, lat05), linewidth=1, color=scol)
      M.plot( (lon05_2, lon2), (lat05, lat2), linewidth=1, color=scol)
      #--------------
    else:
      M.plot( (lon1, lon2), (lat1, lat2), linewidth=1, color=scol)





#-- coastline ---------------
print "coastlines"
M.drawcoastlines(color="k")



figdir  = '/home/utsumi/temp/bams2020'
util.mk_dir(figdir)
figpath = figdir + '/temp.png'
plt.savefig(figpath)
print figpath

sys.exit()
#-- draw cyclone tracks ------
itemp = 1
for iclass in lclass[1:]:
#for iclass in lclass:
  #if iclass < iclass_min: continue
  if (len(dtrack[iclass]) ==  0.0):
    continue
  #-----------
  for track in dtrack[iclass]:
    itemp = itemp + 1
    year = track[0][0]
    mon  = track[0][1]
    day  = track[0][2]
    hour = track[0][3]

    lat1 = track[1][0]
    lon1 = track[1][1]
    lat2 = track[1][2]
    lon2 = track[1][3]

    #---- check region ----------
    if ((lat1 < lllat) or (urlat < lat1)):
      continue
    if ((lon1 < lllon) or (urlon < lon1)):
      continue
    #--------------
    #if iclass ==0:
    #  scol="r"

    if iclass ==1:
      scol="gray"
      #scol="b"
    elif iclass ==2:
      scol="b"
      #scol="r"
    elif iclass ==3:
      scol="r"
#    elif iclass == 4:
#      #scol="gray"
#      scol="r"
 
    #------------------------------------
    if abs(lon1 - lon2) >= 180.0:
      #--------------
      if (lon1 > lon2):
        lon05_1  = 360.0
        lon05_2  = 0.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      elif (lon1 < lon2):
        lon05_1  = 0.0
        lon05_2  = 360.0
        lat05    = lat1 + (lat2 - lat1)/(lon05_1 - lon1 + lon2 - lon05_2)*(lon05_1 - lon1)
      #--------------
      M.plot( (lon1, lon05_1), (lat1, lat05), linewidth=1, color=scol)
      M.plot( (lon05_2, lon2), (lat05, lat2), linewidth=1, color=scol)
      #--------------
    else:
      M.plot( (lon1, lon2), (lat1, lat2), linewidth=1, color=scol)

    #-- text -----------
    if hour in [0,12]:
      xtext, ytext = M(lon1,lat1)
      #plt.text(xtext,ytext-1, "%02d.%02d"%(day,hour) ,fontsize=12, rotation=-90)




#-- coastline ---------------
print "coastlines"
M.drawcoastlines(color="k")

#-- meridians and parallels
parallels = arange(-90.,90,10.)
M.drawparallels(parallels,labels=[1,0,0,0],fontsize=lonlatfontsize)

meridians = arange(0.,360.,10.)
M.drawmeridians(meridians,labels=[0,0,0,1],fontsize=lonlatfontsize,rotation=lonrotation)

#-- title --------------------
stitle  = "%04d/%02d/%02d-%04d/%02d/%02d"\
          %(iDTime.year, iDTime.month, iDTime.day, eDTime.year, eDTime.month, eDTime.day)
axmap.set_title(stitle, fontsize=10.0)

#-- save --------------------
print "save"
#sodir   = "/media/disk2/out/cyclone/exc.track.w.bsttc.JRA55/%s"%(region)
#sodir   = "/home/utsumi/temp"
sodir  = "."
util.mk_dir(sodir)
#soname  = sodir + "/exc.track.w.bsttc.%s.%04d.%02d.%02d-%02d.%02dh.png"%(model, year,mon, iday, eday, thdura)
soname  = sodir + "/test.png"
plt.savefig(soname)
plt.clf()
print soname
plt.clf()
