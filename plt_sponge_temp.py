import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v39/sponge.nc','r')

x = data.variables['x'][:]
y = data.variables['y'][:]
zt = data.variables['zt'][:]
zw = data.variables['zw'][:]

Idamp = data.variables['Idamp'][:,:]
ptemp = data.variables['PTEMP'][:,-1,-1]
eta   = data.variables['ETA'][:,:,:]
salt  = data.variables['SALT'][:,:,:]

# plot - restoring temperature

plt.figure(1,figsize=(5,10))

plt.plot(ptemp,zt,'k')
plt.scatter(ptemp,zt,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(0,12)
plt.ylim(0,4000)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel('PTEMP ($^\circ$C)',fontsize=20)
plt.ylabel('Depth (m)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])


plt.savefig('/short/v45/lxy581/mom6/diag/v39_sponge_temp.png',dpi=600)

plt.show()

