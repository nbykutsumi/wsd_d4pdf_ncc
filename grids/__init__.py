from numpy import *
from bisect import bisect, bisect_left, bisect_right
import numpy as np
import numpy.ma as ma
import sys
import math

def mk_a2gridarea(a1latBnd, a1lonBnd):
    """
    a1latBnd, a1lonBnd in degrees, grid boundaries
    returns in km^2
    based on Oki and Kanae 1997, J. Japan Soc. Hydrol. & Water Resour.
    """
    a1lat0 = a1latBnd[:-1]
    a1lat1 = a1latBnd[1:]
    a1lon0 = a1lonBnd[:-1]
    a1lon1 = a1lonBnd[1:]
    
    Lon0, Lat0 = np.meshgrid(a1lon0, a1lat0)
    Lon1, Lat1 = np.meshgrid(a1lon1, a1lat1)
    
    pi = math.pi
    R = 6378.136
    e = 0.081819191
    e2= 0.006694380
    
    def calc_f(Lat):
        Lat = np.deg2rad(Lat)  # degree --> radian
        return 0.5*np.sin(Lat)/(1-e2*(np.sin(Lat)**2)) + 1/(4*e)*np.log(np.abs((1+e*np.sin(Lat))/(1-e*np.sin(Lat))))  # np.log: natural logarithm
    
    f0 = calc_f(Lat0)
    f1 = calc_f(Lat1)
    
    a2area = pi*R**2*(1-e2)/180 * (f1 - f0) * np.abs(Lon1-Lon0) # km2
    return a2area

def unfold2d( aSrc ):
    '''
    * unfold sub-domains
    from div_domain.py of CoreFrame, by H.Kim @ IIS, Univ. of Tokyo
    '''

    shape       = np.array( aSrc.shape, 'int' )
    outershp    = shape[-2:]
    innershp    = shape[:-2]

    nY          = outershp[0] * innershp[-2]
    nX          = outershp[1] * innershp[-1]

    toshape     = innershp[:-2].tolist() + [nY, nX]
    transAxis   = np.array( [0,2,1,3], 'int' )


    if len( shape ) > 2:

        additionalaxes  = list(range( len(shape)-4))

        transAxis       = ( transAxis + len(additionalaxes) ).tolist()
        idx             = 0

        for i in additionalaxes[::-1]:
            transAxis.insert(idx, i)


    return aSrc.transpose(*transAxis).reshape(*toshape)


def fold2d( aSrc, bndShp ):
    '''
    * divide n-dimensional map into small rectagular domains
    from div_domain.py of CoreFrame, by H.Kim @ IIS, Univ. of Tokyo

    aSrc    : input ndarray
    bndShp  : sub-domain shape

    ex) fold2d( ar(t,180,360), (90,180) ).shape   = (t, 2, 2,   90, 180)
        fold2d( ar(t,180,360), ( 2,  2) ).shape   = (t, 90,180, 2,  2)
    '''

    shape       = np.array( aSrc.shape, 'int' )
    outershp    = np.array( bndShp, 'int' )
    innershp    = shape[-len(outershp):] / outershp

    transAxis   = np.array( [0,2,1,3], 'int' )

    if len( shape ) > 2:

        additionalaxes  = list(range( len(shape) - 2))

        transAxis       = ( transAxis + len(additionalaxes) ).tolist()
        idx             = 0

        for i in additionalaxes[::-1]:
            transAxis.insert(idx, i)

    toshape     = shape[:-len(outershp)].tolist()   \
                + [outershp[0], innershp[0], outershp[1], innershp[1]]


    return aSrc.reshape( *toshape ).transpose( *transAxis )


def karnel_pooling_map2D_global(ain, dy, dx, func=None, miss_in=-9999, miss_out=-9999, cover_poles=True):
    ny,nx = ain.shape
    ndup =  (2*dy+1)*(2*dx+1) # number of duplication

    dy,dx = abs(dy),abs(dx)
  
    if type(ain) is np.ma.core.MaskedArray: 
        ain   = ain.filled(miss_in)

    a3tmp = np.full((ndup,ny,nx),miss_in, dtype=ain.dtype)

    ldyx = [[y,x] for y in range(-dy,dy+1)
                  for x in range(-dx,dx+1)]

    for idup, (idy,idx) in enumerate(ldyx):
        a3tmp[idup] = shift_map_global(ain, idy, idx, cover_poles=cover_poles)

    a3tmp = ma.masked_equal(a3tmp, miss_in)
    print(a3tmp.shape)
    if func =='sum':
        a2out = a3tmp.sum(axis=0)
    elif func=='mean':
        a2out = a3tmp.mean(axis=0)
    elif func=='max':
        a2out = a3tmp.max(axis=0)
    elif func=='min':
        a2out = a3tmp.min(axis=0)
    else:
        print('check func', func)
        sys.exit()

    a2out = a2out.filled(miss_out)
    return a2out


def karnel_pooling_map2D_regional(ain, dy, dx, func=None, miss_in=-9999, miss_out=-9999):
    ny,nx = ain.shape
    ndup =  (2*dy+1)*(2*dx+1) # number of duplication

    dy,dx = abs(dy),abs(dx)
  
    if type(ain) is np.ma.core.MaskedArray: 
        ain   = ain.filled(miss_in)

    a3tmp = np.full((ndup,ny,nx),miss_in, dtype=ain.dtype)

    ldyx = [[y,x] for y in range(-dy,dy+1)
                  for x in range(-dx,dx+1)]

    for idup, (idy,idx) in enumerate(ldyx):
        a3tmp[idup] = shift_map_regional(ain, idy, idx, fill_value=miss_in)

    a3tmp = ma.masked_equal(a3tmp, miss_in)
    if func =='sum':
        a2out = a3tmp.sum(axis=0)
    elif func=='mean':
        a2out = a3tmp.mean(axis=0)
    elif func=='max':
        a2out = a3tmp.max(axis=0)
    elif func=='min':
        a2out = a3tmp.min(axis=0)
    else:
        print('check func', func)
        sys.exit()

    a2out = a2out.filled(miss_out)
    return a2out






def shift_map_global(a2in, dy, dx, cover_poles=True, fill_value=np.nan):
    ny,nx = a2in.shape
    a2large = expand_map_global_2d(a2in, abs(dy), abs(dx), cover_poles=cover_poles, fill_value=fill_value)

    return a2large[abs(dy)-dy:abs(dy)-dy+ny, abs(dx)-dx:abs(dx)-dx+nx]

def shift_map_regional(a2in, dy, dx, fill_value=-9999.):
    ny,nx = a2in.shape
    a2out = np.full([ny,nx], fill_value, dtype=a2in.dtype)
    if dy <0:
        oy0=0; oy1=-abs(dy)
        iy0=abs(dy); iy1=None
    elif dy>0:
        oy0=dy; oy1=None
        iy0=0; iy1=-dy
    else:
        oy0=0; oy1=None
        iy0=0; iy1=None

    if dx <0:
        ox0=0; ox1=-abs(dx)
        ix0=abs(dx); ix1=None
    elif dx >0:
        ox0=dx; ox1=None
        ix0=0; ix1=-dx
    else:
        ox0=0; ox1=None
        ix0=0; ix1=None

    a2out[oy0:oy1,ox0:ox1] = a2in[iy0:iy1,ix0:ix1]
    return a2out


def mk_mask_BBox(a1lat, a1lon, BBox, miss=0):
    """
    Only for global map
    """

    [[lat_min, lon_min],[lat_max, lon_max]] = BBox
    dlon      = (a1lon[1] - a1lon[0])
    dlat      = (a1lat[1] - a1lat[0])
    lon_first = a1lon[0] - 0.5*dlon
    lat_first = a1lat[0] - 0.5*dlat
    lon_last  = a1lon[-1] + 0.5*dlon
    lat_last  = a1lat[-1] + 0.5*dlat

    #--- xmin ----------
    if (lon_first <= lon_min):
        if (lon_min <= lon_last):
            xmin = bisect_right(a1lon+0.5*dlon, lon_min)
        else:
            xmin = bisect_right(a1lon+0.5*dlon, lon_min-lon_last)
    else:
        xmin = bisect_right(a1lon+0.5*dlon, lon_min+lon_last)

    #--- xmax ----------
    if (lon_first <= lon_max):
        if (lon_max <= lon_last):
            xmax = bisect_left(a1lon+0.5*dlon, lon_max)
        else:
            xmax = bisect_left(a1lon+0.5*dlon, lon_max-lon_last)
            
    else:
        xmax = bisect_left(a1lon+0.5*dlon, lon_max+lon_last)

    #--- ymin ----------
    ymin = bisect_right(a1lat+0.5*dlat, lat_min)

    #--- ymax ----------
    ymax = bisect_left(a1lat+0.5*dlat, lat_max)

    ##-----------
    ny = len(a1lat)
    nx = len(a1lon)
    a2regionmask  = ones(nx*ny).reshape(ny, nx)*miss

    if ( xmax < xmin):
        a2regionmask[ymin:ymax+1, xmin: nx] = 1.0
        a2regionmask[ymin:ymax+1, 0:xmax+1] = 1.0
    else:
        a2regionmask[ymin:ymax+1, xmin: xmax+1] = 1.0
    return a2regionmask

def expand_map_global_2d(a2in, dy,dx, cover_poles=True, fill_value=np.nan):
    '''
    expand global map
    input shape (ny,nx)
    output shape: (ny+2dy, nx+2dx)
    '''
    ny,nx = a2in.shape
    dy, dx = abs(dy), abs(dx)     
    #-- Make output array   --
    a2out = np.full([ny+2*dy, nx+2*dx], fill_value).astype(a2in.dtype)
    #-- Fill center box -----
    a2out[dy:dy+ny, dx:dx+nx] = a2in

    #-- Latitude direction --
    if cover_poles==True:
        a2out[:dy,dx:dx+nx]  = np.flipud(np.roll(a2in[1:1+dy,:],  int(0.5*nx), axis=1))
        a2out[-dy:,dx:dx+nx] = np.flipud(np.roll(a2in[-dy-1:-1,:], int(0.5*nx), axis=1))

    else:
        a2out[:dy,dx:dx+nx]  = np.flipud(np.roll(a2in[:dy,:],  int(0.5*nx), axis=1))
        a2out[-dy:,dx:dx+nx] = np.flipud(np.roll(a2in[-dy:,:], int(0.5*nx), axis=1))

    #-- Longitude direction (wrap) --
    if dx !=0:
        a2out[:, :dx]  = a2out[:,-2*dx:-dx]
        a2out[:, -dx:] = a2out[:,dx:2*dx]

    return a2out.astype(a2in.dtype)

def expand_map_regional_2d(a2in, dy,dx, miss_fill):
    '''
    expand regional map
    input shape (ny,nx)
    output shape: (ny+2dy, nx+2dx)
    '''
    ny,nx = a2in.shape
     
    #-- Make output array   --
    a2out = np.full([ny+2*dy, nx+2*dx], miss_fill).astype(a2in.dtype)

    #-- Fill center box -----
    a2out[dy:dy+ny, dx:dx+nx] = a2in

    return a2out.astype(a2in.dtype)


#def shift_map(a2in, dy, dx, miss):
#    ny, nx = a2in.shape
#    a2out  = ones([ny,nx])*miss
#
#    if dy>=0:
#        if dx==0:
#            a2out[dy:,:] = a2in[:ny-dy,:]
#        else:
#            a2out[dy:,dx:] = a2in[:ny-dy,:-dx]
#            a2out[dy:,:dx] = a2in[:ny-dy,-dx:]
#
#    else:
#        if dx==0:
#            a2out[:dy,:] = a2in[-dy:,:]
#        else:
#            a2out[:dy,dx:] = a2in[-dy:,:-dx]
#            a2out[:dy,:dx] = a2in[-dy:,-dx:]
#
#    return a2out

     
