
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

#data = nc.Dataset('/short/v45/lxy581/mom6/work/so_mom6_v4/prog.nc','r')
data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/prog.nc','r')

xq = data.variables['xq'][:]   # u
yh = data.variables['yh'][:]   # u,temp
xh = data.variables['xh'][:]   # v,temp
yq = data.variables['yq'][:]   # v
zl = data.variables['zl'][:]

temp4 = data.variables['temp'][:,:,:,:]

# time-mean
temp3 = np.nanmean(temp4,axis=0)

# zonal-mean
temp2 = np.nanmean(temp3,axis=-1)

# plot

#plt.figure(1,figsize=(5,12))
plt.figure(1,figsize=(8,5))

t_level = np.arange(0,12+0.1,0.1)
t_ticks = np.arange(0,12+4,4)

#plt.subplot(311)
pc3 = plt.contourf(yh,zl,temp2,cmap=plt.cm.jet,levels=t_level,extend='both')
c3 = plt.colorbar(pc3,ticks=t_ticks)
plt.gca().invert_yaxis()
# later versions
plt.gca().set_xticks(np.arange(-1250,1250 + 500,500))
# v2 only
#plt.gca().set_xticks(np.arange(0,2500 + 500,500))
plt.gca().set_xticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_yticks(np.arange(0,4000 + 1000,1000))
plt.gca().set_yticklabels(['0','1','2','3','4'])
#plt.gca().set_position([0.15,0.7,0.65,0.2])
plt.gca().set_position([0.15,0.2,0.65,0.6])
c3.ax.set_position([0.85,0.2,0.03,0.6])
c3.set_label('temp ($^\circ$C)',y=0.5,fontsize=16)
c3.ax.tick_params(labelsize=16)
plt.xlabel('Y (km)',fontsize=16)
plt.ylabel('Depth (km)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

#plt.savefig('/short/v45/lxy581/mom6/diag/v2_ref_merid_5.png',dpi=600)
#plt.savefig('/short/v45/lxy581/mom6/diag/v4_work.png',dpi=600)
plt.savefig('/short/v45/lxy581/mom6/diag/v33_merid_temp_y50.png',dpi=600)

plt.show()





