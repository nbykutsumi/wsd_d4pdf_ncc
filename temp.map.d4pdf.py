# %%
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from   mpl_toolkits.basemap import Basemap
#%matplotlib inline
import d4PDF
from datetime import datetime, timedelta
from numpy import ma
import util

scen = 'HPB'
ens  = 1
dbbaseDir = '/home/utsumi/mnt/lab_work/hk01/d4PDF_GCM'
dbsfc = d4PDF.snp_6hr_2byte(vtype='sfc', dbbaseDir=dbbaseDir)
dbatm = d4PDF.snp_6hr_2byte(vtype='atm', dbbaseDir=dbbaseDir)

#iDTime = datetime(2010,1,1)
#eDTime = datetime(2010,1,19,18)
#lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(hours=6))
DTime = datetime(2010,1,1,18)

vname = 'U850'
a2dat = dbatm.load_6hr(vname=vname, scen=scen, ens=ens, DTime=DTime)
a2dat = ma.masked_equal(a2dat, -9999.)

a1lat = d4PDF.Lat()
a1lon = d4PDF.Lon()
a2lon, a2lat = np.meshgrid(a1lon, a1lat)

#[[lllat,lllon],[urlat,urlon]] = [[0,100],[45,180]]
[[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,360]]
figmap = plt.figure(figsize=(6,4))
axmap  = figmap.add_axes([0.1,0.1,0.8,0.8])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)

im  =  M.pcolormesh( a2lon, a2lat, a2dat)
plt.colorbar(im)
M.drawcoastlines(color="k")
#-- Merdians and Parallels --
M.drawmeridians(np.arange(-180,180+1,15), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50, yoffset=2)
M.drawparallels(np.arange(-60,60+1,15), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')



figpath = '/home/utsumi/temp/bams2020/temp.%s.png'%(vname)
plt.savefig(figpath)
plt.show()
print figpath
# %%
