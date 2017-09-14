import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/prog_rho.nc','r')

xq = data.variables['xq'][:]   # u
yh = data.variables['yh'][:]   # u,temp
xh = data.variables['xh'][:]   # v,temp
yq = data.variables['yq'][:]   # v
rhol = data.variables['02_l'][0:-1]

# ATTENTION: The last number in rhol is very misleading as it it the average of the maximum and minimum density. It has to be excluded otherwise the plots would be very weird. Correspondingly, when reading h_rho and vh_rho, I also exclude the last index for the second dimension (02_l).

h_rho4  = data.variables['h_rho'][:,0:-1,:,:]
vh_rho4 = data.variables['vh_rho'][:,0:-1,:,:]
vh_rho4 = vh_rho4.filled(np.nan)
h_rho4[h_rho4 > 4000.] = 0.0

# time-mean
h_rho3  = np.nanmean(h_rho4,axis=0)
vh_rho3 = np.nanmean(vh_rho4,axis=0)

# zonal-mean
h_rho2  = np.nanmean(h_rho3,axis=-1)
vh_rho2 = np.nanmean(vh_rho3,axis=-1)

# the density larger seem to be assigned the last (largest) layer thickness
[nz,ny] = np.shape(h_rho2)
for j in xrange(0,ny):
    for k in xrange(nz-1,0,-1):
        if h_rho2[k,j] == h_rho2[k-1,j]:
            h_rho2[k,j] = 0.0

h_rho2[h_rho2<2.5]=np.nan
vh_rho2[np.isnan(h_rho2)==True]=np.nan

h_sum   = np.nansum(h_rho2,axis=0)

# plot

plt.figure(1,figsize=(8,5))

h_level = np.arange(0,500+5,5)
h_ticks = np.arange(0,500+100,100)

pc3 = plt.contourf(yh,rhol,h_rho2,cmap=plt.cm.jet,levels=h_level)
#pc3 = plt.contourf(yh,rhol,h_test,cmap=plt.cm.jet,levels=h_level)
c3 = plt.colorbar(pc3,ticks=h_ticks)

plt.gca().set_xticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_xticklabels(['0','500','1000','1500','2000','2500'])
#plt.gca().set_ylim([1035.2,1037.4])
plt.gca().set_ylim([1036.0,1037.2])
plt.gca().set_yticks(np.arange(1036.0,1037.2 + 0.3,0.3))
plt.gca().set_yticklabels(['36.0','36.3','36.6','36.9','37.2'])
plt.gca().invert_yaxis()
plt.gca().set_position([0.15,0.2,0.65,0.6])
c3.ax.set_position([0.85,0.2,0.03,0.6])
c3.set_label('Layer thicknesses (m)',y=0.5,fontsize=16)
c3.ax.tick_params(labelsize=16)
plt.xlabel('Y (km)',fontsize=16)
plt.ylabel('Potential density ($kg\ m^{-3}$)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v33_h_rho_y50.png',dpi=600)

plt.figure(2,figsize=(8,5))

vh_level = np.arange(-0.06e+6,0.06e+6+0.003e+6,0.003e+6)
vh_ticks = np.arange(-0.06e+6,0.06e+6+0.03e+6,0.03e+6)

pc3 = plt.contourf(yq,rhol,vh_rho2,cmap=plt.cm.jet,levels=vh_level,extend='both')
c3 = plt.colorbar(pc3,ticks=vh_ticks)
plt.gca().set_xticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_xticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_ylim([1036.0,1037.2])
plt.gca().set_yticks(np.arange(1036.0,1037.2 + 0.3,0.3))
plt.gca().set_yticklabels(['36.0','36.3','36.6','36.9','37.2'])
plt.gca().invert_yaxis()
plt.gca().set_position([0.15,0.2,0.65,0.6])
c3.ax.set_position([0.85,0.2,0.03,0.6])
c3.ax.set_yticklabels(['-0.06','-0.03','0','0.03','0.06'])
c3.set_label('Meridional transport (Sv)',y=0.5,fontsize=16)
c3.ax.tick_params(labelsize=16)
plt.xlabel('Y (km)',fontsize=16)
plt.ylabel('Potential density ($kg\ m^{-3}$)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v33_vh_rho_y50.png',dpi=600)

plt.show()






