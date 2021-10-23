# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pygrib
import os, sys
%matplotlib inline
import d4PDF
from datetime import datetime, timedelta

#srcdir = '/home/utsumi/mnt/lab_work/hk02/BAMS.EEE.2020/d4PDF_GCM.test/HPB_NAT/m100/201001'
#
##srcpath= srcdir + '/sfc_souseid_avr_day_HPB_NAT_m100_201001.grib'
##grbs = pygrib.open(srcpath)
##grbs.seek(0)
##grb = grbs.select(indicatorOfParameter=255)
#
#srcpath= srcdir + '/sfc_avr_mon_HPB_NAT_m100_201001.grib'
#grbs = pygrib.open(srcpath)
#grbs.seek(0)
#grb = grbs.select(indicatorOfParameter=254)  # 254: TGEF
#data = grb[0].values

db = d4PDF.avr_mon_320x640(vtype='sfc')
data =db.load_ave_mon(vname='TGEF', scen='HPB', ens=1, Year=2010, Mon=1)

print ''
print data
plt.imshow(data,origin='lower')
plt.colorbar()
#data = grbs.select(name='TMPGRD')
#print data

# %%
