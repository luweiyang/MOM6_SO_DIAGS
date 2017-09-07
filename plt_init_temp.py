import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

#data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v3/temp_salt_z.nc','r')
data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v3_1/indata_so_v3.nc','r')

x = data.variables['x'][:]
y = data.variables['y'][:]
zt = data.variables['zt'][:]

itemp = data.variables['itemp'][:,0,0]

# plot - initial temperature

plt.figure(1,figsize=(5,10))

plt.plot(itemp,zt,'k')
plt.scatter(itemp,zt,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(0,12)
plt.ylim(0,4000)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('Initial temp ($^\circ$C)',fontsize=20)
plt.ylabel('Depth (m)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])

plt.savefig('/short/v45/lxy581/mom6/diag/v3_1_init_temp.png',dpi=600)

plt.show()
