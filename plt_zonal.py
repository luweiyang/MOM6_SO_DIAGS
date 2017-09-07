import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/work/so_mom6_v11/prog.nc','r')
#data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v2/output009/prog.nc','r')

xq = data.variables['xq'][:]   # u
yh = data.variables['yh'][:]   # u,temp
xh = data.variables['xh'][:]   # v,temp
yq = data.variables['yq'][:]   # v
zl = data.variables['zl'][:]

#time index
ti = 0

#lat-y index
yi = -10

u = data.variables['u'][ti,:,yi,:]
v = data.variables['v'][ti,:,yi,:]
temp = data.variables['temp'][ti,:,yi,:]

#u = np.transpose(u)

# plot

plt.figure(1,figsize=(5,12))

u_level = np.arange(-1,1+0.1,0.1)
u_ticks = np.arange(-1,1+0.5,0.5)
t_level = np.arange(0,12+0.1,0.1)
t_ticks = np.arange(0,12+4,4)

plt.subplot(311)
pc3 = plt.contourf(xh,zl,temp,cmap=plt.cm.jet,levels=t_level,extend='both')
c3 = plt.colorbar(pc3,ticks=t_ticks)
plt.gca().invert_yaxis()
# later versions
plt.gca().set_xticks(np.arange(-2000,2000 + 800,800))
# v2 only
#plt.gca().set_xticks(np.arange(0,4000 + 800,800))
plt.gca().set_xticklabels(['0','800','1600','2400','3200','4000'])
plt.gca().set_yticks(np.arange(0,4000 + 1000,1000))
plt.gca().set_yticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.1,0.7,0.7,0.2])
c3.ax.set_position([0.85,0.7,0.03,0.2])

plt.subplot(312)
pc1 = plt.contourf(xq,zl,u,cmap=plt.cm.bwr,levels=u_level,extend='both')
c1 = plt.colorbar(pc1,ticks=u_ticks)
plt.gca().invert_yaxis()
plt.gca().set_xticks(np.arange(-2000,2000 + 800,800))
#plt.gca().set_xticks(np.arange(0,4000 + 800,800))
plt.gca().set_xticklabels(['0','800','1600','2400','3200','4000'])
plt.gca().set_yticks(np.arange(0,4000 + 1000,1000))
plt.gca().set_yticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.1,0.4,0.7,0.2])
c1.ax.set_position([0.85,0.4,0.03,0.2])

plt.subplot(313)
pc2 = plt.contourf(xh,zl,v,cmap=plt.cm.bwr,levels=u_level,extend='both')
c2 = plt.colorbar(pc2,ticks=u_ticks)
plt.gca().invert_yaxis()
plt.gca().set_xticks(np.arange(-2000,2000 + 800,800))
#plt.gca().set_xticks(np.arange(0,4000 + 800,800))
plt.gca().set_xticklabels(['0','800','1600','2400','3200','4000'])
plt.gca().set_yticks(np.arange(0,4000 + 1000,1000))
plt.gca().set_yticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.1,0.1,0.7,0.2])
c2.ax.set_position([0.85,0.1,0.03,0.2])

#plt.savefig('/short/v45/lxy581/mom6/diag/v2_ref_zonal_9.png',dpi=600)
plt.savefig('/short/v45/lxy581/mom6/diag/v11_zonal_spg.png',dpi=600)

plt.show()





