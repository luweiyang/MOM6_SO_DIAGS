# TOTAL, BAROTROPIC, BAROCLINIC TRANSPORT ALONG X-AXIS

import numpy as np 
import netCDF4 as nc

from datetime import datetime
start = datetime.now()

# horizontal resolution (km)
dx = 10.0

# bottom index
b5  = nc.Dataset('/short/v45/lxy581/mom6/diag/find_b5.nc','r')
bot = b5.variables['b5low'][:,:]

# example file
examp = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1/output015/prog.nc','r')
xh = examp.variables['xh'][:]
yh = examp.variables['yh'][:]
zl = examp.variables['zl'][:]

nx = np.size(xh)
ny = np.size(yh)
nz = np.size(zl)

# vertical coordinate file
coord = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1/output015/Vertical_coordinate.nc','r')

# layer thicknesses in m
ds = coord.variables['ds'][:]

fac_bottom = np.zeros((nz,ny,nx))
fac_water  = np.zeros((nz,ny,nx))

ds_xyz = np.zeros((nz,ny,nx))
ds_bottom = np.zeros((nz,ny,nx))
ds_water = np.zeros((nz,ny,nx))

for j in xrange(ny):
  for i in xrange(nx):
    fac_bottom[bot[j,i]-1,j,i]  = 1.0
    fac_water[0:bot[j,i],j,i] = 1.0
    ds_xyz[:,j,i] = ds

ds_water  = ds_xyz * fac_water         # 3D xyz
ds_bottom = np.nansum(ds_water,axis=0) # total layer thickness from bottom to top, 2D xy 

# TRANSPORT, per year, along x-axis
tt_yr_x = np.zeros((1,nx))
bt_yr_x = np.zeros((1,nx))
bc_yr_x = np.zeros((1,nx))

tt_yr_xt = np.full((1,nx),np.nan)
bt_yr_xt = np.full((1,nx),np.nan)
bc_yr_xt = np.full((1,nx),np.nan)

# TRANSPORT, per unit time, along x-axis
tt_x  = np.zeros((1,nx))
tt_xt = np.full((1,nx),np.nan)
bt_x  = np.zeros((1,nx))
bt_xt = np.full((1,nx),np.nan)
bc_x  = np.zeros((1,nx))
bc_xt = np.full((1,nx),np.nan)

# LAST 10 YEARS
for opt in xrange(14,15+1,1):
  print '--------'
  print opt
  print '/g/data1/v45/lxy581/payu/archive/exp1_w1/output0%02d/' % opt
  d_opt = '/g/data1/v45/lxy581/payu/archive/exp1_w1/output0%02d/prog.nc' % opt

  data  = nc.Dataset(d_opt,'r')

  t  = data.variables['Time'][:]
  nt = np.size(t)

  for k in xrange(0,nt):
    # 3D zonal velocity over the domain
    u_xyz = data.variables['u'][k,:,:,:] * fac_water  
    # find bottom velocity
    u_bot_3D = u_xyz * fac_bottom
    u_bot = np.nansum(u_bot_3D,axis=0)  # 2D xy
    # find ocean velocity (above bathymetry)
    u_water = u_xyz * fac_water      # 3D xyz

    # bottom u * h
    ubt_xy = u_bot * ds_bottom * (dx * 1.e+3)
    # barotropic transport along x-axis
    ubt_x = np.nansum(ubt_xy,axis=0)

    # water u * h
    utt_xyz = u_water * ds_water * (dx * 1.e+3)
    # total transport along x-axis
    utt_xy = np.nansum(utt_xyz,axis=0) # sum along z
    utt_x = np.nansum(utt_xy,axis=0)  # sum along y
      
    # baroclinic transport along x-axis
    ubc_xy = utt_xy - ubt_xy
    ubc_x = np.nansum(ubc_xy,axis=0) 

    # append
    tt_xt[0,:] = utt_x
    bt_xt[0,:] = ubt_x
    bc_xt[0,:] = ubc_x
    tt_x = np.append(tt_x,tt_xt,0)   
    bt_x = np.append(bt_x,bt_xt,0)   
    bc_x = np.append(bc_x,bc_xt,0)   

tt_x = tt_x[1:,:] 
bt_x = bt_x[1:,:] 
bc_x = bc_x[1:,:]

# count years
yr = np.shape(tt_x)[0]/73
print yr

# calculate annual transport 
for i in xrange(0,yr):
  tt_x_annual = np.nanmean(tt_x[73*i:73*(i+1),:],axis=0)
  bt_x_annual = np.nanmean(bt_x[73*i:73*(i+1),:],axis=0)
  bc_x_annual = np.nanmean(bc_x[73*i:73*(i+1),:],axis=0)

  tt_yr_xt[0,:] = tt_x_annual 
  bt_yr_xt[0,:] = bt_x_annual 
  bc_yr_xt[0,:] = bc_x_annual 

  tt_yr_x = np.append(tt_yr_x,tt_yr_xt,0) 
  bt_yr_x = np.append(bt_yr_x,bt_yr_xt,0) 
  bc_yr_x = np.append(bc_yr_x,bc_yr_xt,0) 

tt_yr_x = tt_yr_x[1:,:]
bt_yr_x = bt_yr_x[1:,:]
bc_yr_x = bc_yr_x[1:,:]

time_diff = datetime.now()-start

# save
fileobj = nc.Dataset('/short/v45/lxy581/mom6/diag/domain_transport_10yr_11_add_spg.nc','w')

fileobj.createDimension('year',yr)
fileobj.createDimension('X',nx)
 
yr_var = fileobj.createVariable('year','f',('year',))
x_var = fileobj.createVariable('x','f',('X',))

tt_var = fileobj.createVariable('tt','f',('year','X'),fill_value=-1e+20)
bt_var = fileobj.createVariable('bt','f',('year','X'),fill_value=-1e+20)
bc_var = fileobj.createVariable('bc','f',('year','X'),fill_value=-1e+20)

yr_var[:] = np.arange(yr)
x_var[:] = xh
tt_var[:,:] = tt_yr_x[:,:] # 10 years of zonal transport
bt_var[:,:] = bt_yr_x[:,:] # 10 years of barotropic transport
bc_var[:,:] = bc_yr_x[:,:] # 10 years of baroclinic transport

print 'Finished'

