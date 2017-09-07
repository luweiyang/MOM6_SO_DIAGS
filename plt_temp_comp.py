# compare sponge temp and real nothernmost temp

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v32/sponge.nc','r')
outp = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v32/output010/prog.nc','r')

x = data.variables['x'][:]
y = data.variables['y'][:]
zt = data.variables['zt'][:]
zw = data.variables['zw'][:]

Idamp = data.variables['Idamp'][:,:]
ptemp = data.variables['PTEMP'][:,-1,-1]
eta   = data.variables['ETA'][:,:,:]
salt  = data.variables['SALT'][:,:,:]

temp3 = outp.variables['temp'][:,:,-1,:]
# time-mean
temp2 = np.nanmean(temp3,axis=0)
# zonal-mean
temp = np.nanmean(temp2,axis=-1)

# plot - restoring temperature

plt.figure(1,figsize=(5,10))

plt.plot(ptemp,zt,'k',label='restoring')
plt.scatter(ptemp,zt,color='k')
plt.plot(temp,zt,'r',label='real')
plt.scatter(temp,zt,color='r')
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
plt.legend(loc=[0.45,0])

plt.savefig('/short/v45/lxy581/mom6/diag/v32_temp_comp.png',dpi=600)

plt.show()

