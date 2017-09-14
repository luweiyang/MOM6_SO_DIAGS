# diagnose bottom kinetic energy
# last year - every 5 days' outputs

print 'Running'

from datetime import datetime
start = datetime.now()


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io

data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/prog.nc','r')
b5   = nc.Dataset('/short/v45/lxy581/mom6/diag/find_b5.nc','r')

xh = data.variables['xh'][:]
yh = data.variables['yh'][:]
zl = data.variables['zl'][:]
t  = data.variables['Time'][:]

nz = np.size(zl)
nx = np.size(xh)
ny = np.size(yh)
nt = np.size(t)

b5up = b5.variables['b5up'][:,:]
b5lo = b5.variables['b5low'][:,:]

# lowest 500m-averaged velocity
u5 = np.zeros((1,ny,nx))
v5 = np.zeros((1,ny,nx))
u5t = np.full((1,ny,nx),np.nan)
v5t = np.full((1,ny,nx),np.nan)
fac = np.full((nz,ny,nx),np.nan)

for j in xrange(ny):
  for i in xrange(nx):
    fac[b5up[j,i]:b5lo[j,i]+1,j,i] = 1.0

for k in xrange(0,nt):
  print k
  u = data.variables['u'][k,:,:,:]
  v = data.variables['v'][k,:,:,:]

  uu = u*fac
  vv = v*fac

  u5t[0,:,:] = np.nanmean(uu,axis=0)
  v5t[0,:,:] = np.nanmean(vv,axis=0)
 
  u5 = np.append(u5,u5t,0)
  v5 = np.append(v5,v5t,0)   

u5 = u5[1:,:,:]
v5 = v5[1:,:,:]
nt = np.size(u5,0)

u5m = np.nanmean(u5,axis=0)
v5m = np.nanmean(v5,axis=0)

u5p = u5 - u5m
v5p = v5 - v5m

MKEb5 = 0.5*(u5m**2 + v5m**2)
TKEb5 = 0.5*(u5**2 + v5**2)
mTKEb5 = np.nanmean(TKEb5,axis=0)
mTKEb50 = np.log10(mTKEb5)

# bottom EKE
EKEb5  = 0.5*(u5p**2 + v5p**2)
mEKEb5 = np.nanmean(EKEb5,axis=0)

# bottom speed
spdb5  = np.sqrt(u5**2 + v5**2)
mspdb5 = np.nanmean(spdb5,axis=0)

'''
# save
fileobj = nc.Dataset('/short/v45/lxy581/mom6/diag/bot_500_SO_10_EKE_spd_170823.nc','w')

fileobj.createDimension('yh',len(yh))
fileobj.createDimension('xh',len(xh)) 
fileobj.createDimension('time',nt)
 
yh_var = fileobj.createVariable('yh','f',('yh',))
xh_var = fileobj.createVariable('xh','f',('xh',))
time_var = fileobj.createVariable('time','f',('time'),)

u5_var = fileobj.createVariable('u5','f',('time','yh','xh'),fill_value=-1e+20)
v5_var = fileobj.createVariable('v5','f',('time','yh','xh'),fill_value=-1e+20)
u5p_var = fileobj.createVariable('u5p','f',('time','yh','xh'),fill_value=-1e+20)
v5p_var = fileobj.createVariable('v5p','f',('time','yh','xh'),fill_value=-1e+20)
u5m_var = fileobj.createVariable('u5m','f',('yh','xh'),fill_value=-1e+20)
v5m_var = fileobj.createVariable('v5m','f',('yh','xh'),fill_value=-1e+20)

mEKEb5_var = fileobj.createVariable('mEKEb5','f',('yh','xh'),fill_value=-1e+20)
mspdb5_var = fileobj.createVariable('mspdb5','f',('yh','xh'),fill_value=-1e+20)

yh_var[:] = yh[:]
xh_var[:] = xh[:]
time_var[:] = t[:]

u5_var[:,:,:] = u5[:,:,:]
v5_var[:,:,:] = v5[:,:,:]
u5p_var[:,:,:] = u5p[:,:,:]
v5p_var[:,:,:] = v5p[:,:,:]
u5m_var[:,:] = u5m[:,:]
v5m_var[:,:] = v5m[:,:]

mEKEb5_var[:,:] = mEKEb5[:,:]
mspdb5_var[:,:] = mspdb5[:,:]
'''
plt.figure(1,figsize=(8,5))

# log - nikurashin and ferrari 2011
kelev = np.arange(-4,-2+0.05,0.05)
ketik = np.arange(-4,-2+1,1)

X,Y=np.meshgrid(xh,yh)

pc = plt.contourf(X,Y,mTKEb50,cmap=plt.cm.jet,levels=kelev,extend='both')
cb = plt.colorbar(pc,ticks=ketik)
plt.gca().set_yticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_yticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_xticks(np.arange(-2000,2000 + 1000,1000))
plt.gca().set_xticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.15,0.2,0.65,0.6])
cb.ax.set_position([0.85,0.2,0.03,0.6])
cb.set_label('Bottom KE ($log_{10}(m^{2}/s^{2})$)',y=0.5,fontsize=16)
cb.ax.tick_params(labelsize=16)
plt.xlabel('X (km)',fontsize=16)
plt.ylabel('Y (km)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v33_b5_mTKE_y50.png',dpi=600)

plt.show()

time_diff = datetime.now()-start

print 'Finished'

