import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v2/output009/prog.nc','r')

x = data.variables['xh'][:]
y = data.variables['yh'][:]
zl = data.variables['zl'][:]
t = data.variables['Time'][:]

temp0 = data.variables['temp'][-1,:,-1,:]
temp = np.nanmean(temp0,axis=1) 

# plot - temperature, northern boundary

plt.figure(1,figsize=(5,10))

plt.plot(temp,zl,'k')
plt.scatter(temp,zl,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(0,12)
plt.ylim(0,4000)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('TEMP ($^\circ$C)',fontsize=20)
plt.ylabel('Depth (m)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])

plt.savefig('/short/v45/lxy581/mom6/diag/v2_north_temp.png',dpi=600)

plt.show()
