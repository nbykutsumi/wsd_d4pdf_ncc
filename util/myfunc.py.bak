from numpy import *
#********************************
def hjoin_str(lstr):
  nline = len(lstr[0].strip().split("\n"))
  ddat = {}
  for i in range(nline):
    ddat[i] = []
  #---------------------
  for strdat in lstr:
    lines  = strdat.strip().split("\n")
    for i, line in enumerate(lines):
      ddat[i].append(line)
  #---------------------
  lout = []
  for i in range(nline):
    lout.append(",".join(ddat[i]))
  #---------------------
  sout = "\n".join(lout)
  return sout 
#********************************
def a2_to_csv(a2in):
  lline = [",".join( map(str,line)) for line in a2in]
  sout  = "\n".join(lline).strip()
  return sout
#********************************
def divide_a1latlon(a1lat_in, a1lon_in, ndiv):
  nlat_in = len(a1lat_in)
  nlon_in = len(a1lon_in)
  #*******************:
  # make a1lon_out
  #---------------------
  nlon_out  = nlon_in * ndiv
  a1lon_out = ones([nlon_out], float32)*-9999.
  #---------------------
  # find a1lon_out(1)
  # back to west
  #----------------------
  dlon_in   = a1lon_in[1] - a1lon_in[0]
  dlon_out  = float(dlon_in / ndiv)
  for i in range(ndiv):
    lon = 0.5*(a1lon_in[0] + a1lon_in[1]) - 0.5*dlon_out - i*dlon_out 
    if (lon <0.):
      a1lon_out[0] = lon + dlon_out
      break
    if i == ndiv:
      a1lon_out[0] = lon
  #--- lon: i>2 ---
  for i in range(1,nlon_out):
    a1lon_out[i] = a1lon_out[0] + i*dlon_out
  #*******************:
  # make a1lat_out
  #----------------------
  if (a1lat_in[0] ==-90.0)&(a1lat_in[-1]==90.0):
    nlat_out  = nlat_in * ndiv
    a1lat_out = ones([nlat_out],float32)
    #-----------------
    # for early lats
    #-----------------
    inbound_n = 0.5*(a1lat_in[1] + ( -90.0) )
    dlat_out  = (inbound_n - (-90.0)) / (ndiv -0.5)
    a1lat_out[0] = -90.0
    for i in range(1,ndiv):
      a1lat_out[i] = -90.0 + i*dlat_out
    #-----------------
    # for ordinary grids
    #-----------------
    for iin in range(1,nlat_in-1):
      inbound_s = 0.5*( a1lat_in[iin]   +  a1lat_in[iin-1] )
      inbound_n = 0.5*( a1lat_in[iin+1] +  a1lat_in[iin]   )
      dlat_out  = (inbound_n - inbound_s) / ndiv
      for i in range(ndiv):
        a1lat_out[ iin *ndiv + i] = inbound_s + 0.5*dlat_out + i*dlat_out
    #-----------------
    # for early lats
    #-----------------
    inbound_s = 0.5*(a1lat_in[-1] + a1lat_in[-2])
    dlat_out  = (90.0 - inbound_s) / (ndiv -0.5)
    for i in range(0,ndiv):
      a1lat_out[nlat_out-1 -i] = 90.0 - i*dlat_out
  #----------
  else:
    nlat_out  = nlat_in * ndiv
    a1lat_out = ones([nlat_out],float32)
    #-----------------
    # for early lats
    #-----------------
    inbound_s = -90.0
    inbound_n = 0.5*(a1lat_in[1] + a1lat_in[0])
    dlat_out =  (inbound_n - inbound_s) / ndiv
    for i in range(ndiv):
      a1lat_out[i] = -90.0 + 0.5*dlat_out + i*dlat_out
    #-----------------
    # for ordinary grids
    #-----------------
    for iin in range(1,nlat_in-1):
      print iin
      inbound_s = 0.5*( a1lat_in[iin]   +  a1lat_in[iin-1] )
      inbound_n = 0.5*( a1lat_in[iin+1] +  a1lat_in[iin]   )
      dlat_out  = (inbound_n - inbound_s) / ndiv
      for i in range(ndiv):
        a1lat_out[ iin *ndiv + i] = inbound_s + 0.5*dlat_out + i*dlat_out
    #-----------------
    # for final lats
    #-----------------
    inbound_s = 0.5*(a1lat_in[nlat_in-1] + a1lat_in[nlat_in-1 -1])
    inbound_n = 90.0
    dlat_out  = (inbound_n - inbound_s) / ndiv
    for i in range(ndiv):
        a1lat_out[ (nlat_in -1)*ndiv + i] = inbound_s + 0.5*dlat_out + i*dlat_out
  #-------------
  return a1lat_out, a1lon_out






