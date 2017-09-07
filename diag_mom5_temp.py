# MOM5 offline diag

# 40S temperature vertical profile (zonal mean)
# plot temp(z)

# 40S - 65S surface temperature (zonal mean)
# plot temp(y)

from datetime import datetime
start = datetime.now()

import numpy as np 
import netCDF4 as nc
import scipy.io
import matplotlib.pyplot as plt

grid = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/ocean_grid.nc','r') 
opt  = nc.Dataset('/g/data1/v45/mom/mom01v3/output065/temp_snap.nc','r') 

x = grid.variables['xu_ocean'][:]
y = grid.variables['yu_ocean'][:]
z = opt.variables['st_ocean'][:]

# find the latitude of the Southern Ocean --- 40s~65S
idxy = np.where(np.logical_and(y>=-65,y<=-40))
lat = y[idxy]                         

# find the longitude of the Southern Ocean --- 180W~180E 
idxx = np.where(np.logical_and(x>=-280,x<=80))
lon = x[idxx]

ny = np.size(idxy[0]) # ny = 426
nx = np.size(idxx[0]) # nx = 3600
nz = np.size(z)

# surface
ts = np.zeros((1,ny,nx))
# 40S
t4 = np.zeros((1,nz,nx))

tst = np.full((1,ny,nx),np.nan)
t4t = np.full((1,nz,nx),np.nan)

for opt in xrange(80,92,1):
  print '--------'
  print opt
  dir_opt = '/g/data1/v45/mom/mom01v3/output0%d/temp_snap.nc' % opt

  data = nc.Dataset(dir_opt,'r')
  
  t = data.variables['time'][:]
  print np.shape(t)

  for k in xrange(0,np.size(t),20):
# for k in xrange(3):
    print k
    temp_sf = data.variables['temp'][k,0,idxy[0][:],idxx[0][:]]
    temp_40 = data.variables['temp'][k,:,idxy[0][-1],idxx[0][:]]

    temp_sf = temp_sf.filled(np.nan)
    temp_40 = temp_40.filled(np.nan)

    tst[0,:,:] = temp_sf[:,:]
    t4t[0,:,:] = temp_40[:,:]

    ts = np.append(ts,tst,0) 
    t4 = np.append(t4,t4t,0) 

ts = ts[1:,:,:]
t4 = t4[1:,:,:]

tsm = np.nanmean(ts,axis=0)
t4m = np.nanmean(t4,axis=0)

tsmy = np.nanmean(tsm,axis=1)
t4mz = np.nanmean(t4m,axis=1)

# plot tsmy(y)

plt.figure(1,figsize=(10,6))

plt.plot(tsmy,lat,linewidth=2.0,color='k')
plt.xlim(-3,15)
plt.ylim(-65,-40)
plt.xticks(np.arange(-3,15+3,3))
plt.gca().set_xticklabels(['-3','0','3','6','9','12','15'])
plt.yticks(np.arange(-60,-45+5,5))
plt.gca().set_yticklabels(['60','55','50','45'])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
plt.xlabel('PTEMP ($^\circ$C)',fontsize=20)
plt.ylabel('LAT ($^\circ$S)',fontsize=20)
plt.gca().set_position([0.2,0.2,0.6,0.6])

plt.savefig('/short/v45/lxy581/mom6/diag/diag_mom5_surf_temp.png',dpi=600)

# plot t4mz(z)

plt.figure(2,figsize=(5,10))

plt.plot(t4mz,z,'k')
plt.scatter(t4mz,z,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(0,16)
plt.ylim(0,5000)
plt.xticks(np.arange(4,16+4,4))
plt.gca().set_xticklabels(['4','8','12','16'])
plt.yticks(np.arange(0,5000+1000,1000))
plt.gca().set_yticklabels(['0','1','2','3','4','5'])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('PTEMP ($^\circ$C)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])

plt.savefig('/short/v45/lxy581/mom6/diag/diag_mom5_40S_temp.png',dpi=600)

plt.show()

time_diff = datetime.now()-start
print 'Finished! ^.^'
