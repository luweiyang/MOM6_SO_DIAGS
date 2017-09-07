import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io

# read u
data  = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/prog.nc','r')
# read layer thicknesses
coord = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/Vertical_coordinate.nc','r')
# horizontal resolution in km
dx = 10.

yh = data.variables['yh'][:]
xq = data.variables['xq'][:]
zl = data.variables['zl'][:]
t  = data.variables['Time'][:]

# layer thicknesses in m
ds = coord.variables['ds'][:]

nz = np.size(zl)
nx = np.size(xq)
ny = np.size(yh)
nt = np.size(t)

# zonal- and time-averaged transport - ACC transport
acc  = np.zeros((1,nz))
acct = np.full((1,nz),np.nan)

for k in xrange(0,nt):
  print k

  u  = data.variables['u'][k,:,:,:]  # 3D zonal velocity
  uz = np.nanmean(u,axis=2)          # zonal-mean velocity uz(z,y)
  uy = np.nansum(uz,axis=1)          # meridional sum velocity uy(z)
  uh = uy*ds*dx*1.e+3                # zonal transport uh(z) in m3/s

  acct[0,:] = uh
 
  acc = np.append(acc,acct,0)

acc  = acc[1:,:]
accm = np.nanmean(acc,axis=0) # accm(z)

acc_sum = np.nansum(accm)     # total transport

# acc_sum = 162696873.73828626 m3/s = 162.70 SV

# plot - vertical distribution of zonal transport

plt.figure(1,figsize=(5,10))

plt.plot(accm,zl,'k',linewidth=2)
plt.scatter(accm,zl,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(-2e+6,8e+6)
plt.ylim(0,4000)
plt.xticks(np.arange(-2e+6,8e+6+2e+6,2e+6))
plt.gca().set_xticklabels(['-2','0','2','4','6','8'])
plt.yticks(np.arange(0,4000+1000,1000))
plt.gca().set_yticklabels(['0','1','2','3','4'])
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('Transport (Sv)',fontsize=20)
plt.ylabel('Depth (km)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])

#plt.savefig('/short/v45/lxy581/mom6/diag/v32_y30_acc_170831.png',dpi=600)
plt.savefig('/short/v45/lxy581/mom6/diag/v33_y50_acc_z_170907.png',dpi=600)

plt.show()

