import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/prog.nc','r')

xq   = data.variables['xq'][:]   # u
yh   = data.variables['yh'][:]   # u,temp,rho
xh   = data.variables['xh'][:]   # v,temp,rho
zl   = data.variables['zl'][:]   # v
t    = data.variables['Time'][:]

nz = np.size(zl)
ny = np.size(yh)
nt = np.size(t)

rho3 = np.zeros((1,nz,ny))
rhox = np.full((1,nz,ny),np.nan)

for k in xrange(0,nt):
    print k
    rho  = data.variables['rhoinsitu'][k,:,:,:]

    rhox[0,:,:] = np.nanmean(rho,axis=-1) # zonal-mean
 
    rho3 = np.append(rho3,rhox,0)

rho3 = rho3[1:,:,:]
rho2 = np.nanmean(rho3,axis=0)
prho = rho2 - 1.e+3

plt.figure(1,figsize=(8,5))

rholev = np.arange(35.6,36.8+0.3,0.3)

pc = plt.contour(yh,zl,prho,cmap=plt.cm.jet,levels=rholev)
manual_locations = [(750,300),(400,400),(750,900),(750,1400),(750,2400)]
#plt.clabel(pc,pc.levels[::2],inline=1,fontsize=10,manual=manual_locations)
plt.clabel(pc,inline=1,fontsize=10,manual=manual_locations)
plt.gca().invert_yaxis()
plt.gca().set_xticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_xticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_yticks(np.arange(0,4000 + 1000,1000))
plt.gca().set_yticklabels(['0','1','2','3','4'])
plt.gca().set_position([0.15,0.2,0.70,0.6])
plt.title('In situ density ($kg\ m^{-3}$)',fontsize=16)
plt.xlabel('Y (km)',fontsize=16)
plt.ylabel('Depth (km)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v33_merid_rho_y50.png',dpi=600)

plt.show()

