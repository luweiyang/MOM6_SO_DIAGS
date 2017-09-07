import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v6/topog.nc','r')

x = data.variables['x'][:]
y = data.variables['y'][:]
d = data.variables['depth'][:,:]

# plot - topography

plt.figure(1,figsize=(8,5))

dlevel = np.arange(0,4000+5,5)
dticks = np.arange(0,4000+1000,1000)

pc = plt.contourf(x,y,d,cmap=plt.cm.jet,levels=dlevel) 
cb = plt.colorbar(pc,ticks=dticks)
plt.xticks(np.arange(-2000.e+3,2000.e+3 + 1000.e+3,1000.e+3))
plt.gca().set_xticklabels(['0','1000','2000','3000','4000'])
plt.xlabel('X',fontsize=18)
plt.yticks(np.arange(-1250.e+3,1250.e+3 + 500.e+3,500.e+3))
plt.gca().set_yticklabels(['0','500','1000','1500','2000','2500'])
plt.ylabel('Y',fontsize=18)
plt.gca().tick_params(axis='x',labelsize=14)
plt.gca().tick_params(axis='y',labelsize=14)
plt.gca().set_position([0.15,0.2,0.65,0.6])
cb.ax.set_position([0.85,0.2,0.08,0.6])
cb.ax.tick_params(labelsize=14)
cb.set_label('Depth (m)',y=0.5,fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v6_topog.png',dpi=600)

plt.show()

