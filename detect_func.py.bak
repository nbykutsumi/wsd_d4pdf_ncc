from numpy import *
import os
import grids
import numpy as np
####################################################
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#---------------------------------------------------
def solve_time(stime):
  year = int( stime/10**6 )
  mon  = int( (stime - year*10**6)/10**4 )
  day  = int( (stime - year*10**6 - mon*10**4)/10**2)
  hour = int( (stime - year*10**6 - mon*10**4 - day*10**2) )
  return year, mon, day, hour
#---------------------------------------------------
def fortpos2pyxy(number, nx, miss_int):
  if (number == miss_int):
    iy_py = miss_int
    ix_py = miss_int
  else:
    iy_py = int((number-1.0)/nx)      # iy_py = 0,1,2,..
    ix_py = number - nx*iy_py-1   # ix_py = 0,1,2,..
  #----
  return ix_py, iy_py
#---------------------------------------------------
def box_filtering(a2in, a1y, a1x, dy, dx, func=None, cover_poles=True, miss_in=-9999., miss_out=-9999.):

  if type(a2in) is np.ma.core.MaskedArray:
    a2in = a2in.filled(miss_in)

  a2large = grids.expand_map_global_2d(a2in, dy, dx, cover_poles=cover_poles)
  a1y_large = a1y+dy
  a1x_large = a1x+dx
  ldyx = [[idy,idx] for idy in range(-dy,dy+1) for idx in range(-dx, dx+1)]
  nybox = 2*dy+1
  nxbox = 2*dx+1

  a2container = np.empty([len(a1y), nybox*nxbox], a2in.dtype)
  for i,(idy,idx) in enumerate(ldyx):
    a2container[:, i] = a2large[a1y_large+idy, a1x_large+idx]


  #print a2container


  if func=='sum':
    a1out = ma.masked_equal(a2container, miss_in).sum(axis=1)
  elif func=='mean':
    a1out = ma.masked_equal(a2container, miss_in).mean(axis=1)
  elif func=='min':
    a1out = ma.masked_equal(a2container, miss_in).min(axis=1)
  elif func=='max':
    a1out = ma.masked_equal(a2container, miss_in).max(axis=1)
  else:
    print 'check func',func
    sys.exit()

  return a1out.filled(miss_out)

