# %%
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from   mpl_toolkits.basemap import Basemap
from   datetime import datetime, timedelta
from   numpy import ma
import util
import sys, os
import d4PDF

vname = 'idate'
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/bams2020/XX-HPB_NAT-100/6hr/%s/2010/01'%(vname)
idtime = datetime(2010,1,1,0)
edtime = datetime(2010,1,31,18)
ldtime = util.ret_lDTime(idtime, edtime, timedelta(hours=6))
ny,nx = 320,640

a1lat = d4PDF.Lat()
a1lon = d4PDF.Lon()
[[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,360]]

a2out = np.zeros([ny,nx], 'int32')
a2zero= np.zeros([ny,nx], 'int32')
for dtime in ldtime:
    year,mon,day,hour = dtime.timetuple()[:4]
    srcpath = srcdir + '/%s.%04d%02d%02d%02d.320x640'%(vname,year,mon,day,hour)
    #a2in = np.fromfile(srcpath,'float32').reshape(ny,nx) 
    a2in = np.fromfile(srcpath,'int32').reshape(ny,nx) 
    a2out = a2out + ma.masked_where(a2in>0, a2zero).filled(1).astype('int32')


figmap   = plt.figure()
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
im = M.imshow(ma.masked_equal(a2out,0), origin='lower',vmax=5)
plt.colorbar(im)
M.drawcoastlines()
plt.savefig('/home/utsumi/temp/bams2020/temp.centers.%s.png'%(vname))