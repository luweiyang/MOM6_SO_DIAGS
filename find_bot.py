# bottom 500m index range

from datetime import datetime
start = datetime.now()

import numpy as np 
import netCDF4 as nc

# read topography 
topog = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v6/topog.nc','r')

x = topog.variables['x'][:]
y = topog.variables['y'][:]
d = topog.variables['depth'][:,:]

# read vertical coordinate from outputs
data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v6/output002/prog.nc','r')

xq = data.variables['xq'][:]
yh = data.variables['yh'][:]
zl = data.variables['zl'][:]
xh = data.variables['xh'][:]
yq = data.variables['yq'][:]

nz = np.size(zl)
nx = np.size(xh)
ny = np.size(yh)

b5_low = np.full((ny,nx),np.nan)
b5_up  = np.full((ny,nx),np.nan)

# Note that x and y are 1000 times larger than xh and yh

# To find the layers within 500m above topography for each xh and yh

dpr = 500.0

for j in xrange(ny):
    print j
    for i in xrange(nx):
        for k in np.arange(1,nz):
            if zl[k] > d[j,i] and zl[k-1] < d[j,i]:
               b5_low[j,i] = k-1
            elif zl[nz-1] < d[j,i]:
               b5_low[j,i] = 71
        for m in np.arange(1,nz):
            if zl[m] > d[j,i]-dpr and zl[m-1] < d[j,i]-dpr:
               b5_up[j,i] = m
            elif d[j,i] < dpr:
               b5_up[j,i] = 0

# write
fileobj = nc.Dataset('/short/v45/lxy581/mom6/diag/find_b5.nc','w')
fileobj.createDimension('yh',ny) 
fileobj.createDimension('xh',nx) 
yh_var = fileobj.createVariable('yh','i',('yh',))
xh_var = fileobj.createVariable('xh','i',('xh',))
b5up_var = fileobj.createVariable('b5up','i',('yh','xh'))
b5low_var = fileobj.createVariable('b5low','i',('yh','xh'))
yh_var[:] = yh[:]
xh_var[:] = xh[:]
b5up_var[:,:] = b5_up[:,:]
b5low_var[:,:] = b5_low[:,:]
 
time_diff = datetime.now()-start

print 'Finished'

   





