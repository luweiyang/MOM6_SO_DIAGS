import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
from mpl_toolkits.basemap import Basemap

'''
data = nc.Dataset('/short/v45/lxy581/mom6/input/so_mom6_v2/indata_so_v2.nc','r')

x = data.variables['x'][:]
y = data.variables['y'][:]

sst = data.variables['sst'][:,0]
taux = data.variables['taux'][:,0]

######
# FIGURE 1 - plot surface forcing
fig1 = plt.figure(1,figsize=(10,4))
# wind forcing
ax1 = fig1.add_subplot(111)
ax1.plot(taux,y,linewidth=2.0,color='b')
ax1.set_yticks(np.arange(-1250e+3,1250e+3 + 500e+3,500e+3))
ax1.set_yticklabels(['0','500','1000','1500','2000','2500'])
ax1.set_ylabel('Y (km)',fontsize=20)
ax1.set_xlim(0,0.16)
ax1.set_xticks(np.arange(0,0.16+0.04,0.04))
ax1.set_xticklabels(['0','0.04','0.08','0.12','0.16'])
ax1.xaxis.set_ticks_position('top')
ax1.set_xlabel('Wind forcing',fontsize=20)
ax1.xaxis.set_label_position('top')
ax1.spines['bottom'].set_color('red')
ax1.spines['top'].set_color('blue')
ax1.tick_params(axis='x',colors='blue',labelsize=16)
ax1.tick_params(axis='y',labelsize=16)
ax1.grid(True)
ax1.set_position([0.2,0.2,0.6,0.6])

# surface restoring temperature
ax2 = fig1.add_subplot(111,sharey=ax1,frameon=False)
ax2.plot(sst,y,linewidth=2.0,color='r')
ax2.set_xlim(-4,12)
ax2.set_xticks(np.arange(-4,12+4,4))
ax2.set_xticklabels(['-4','0','4','8','12'])
ax2.xaxis.set_ticks_position('bottom')
ax2.set_xlabel('Surface restoring temp ($^\circ$C)',fontsize=20)
ax2.xaxis.set_label_position('bottom')
ax2.tick_params(axis='x',colors='red',labelsize=16)
ax2.tick_params(axis='y',labelsize=16)
ax2.set_position([0.2,0.2,0.6,0.6])

plt.savefig('/short/v45/lxy581/mom6/diag/v2_forcing.png',dpi=600)
'''

######
# FIGURE 2 - COMPARE 
# (1) surface restoring; with 
# (2) surface zonal- and time- mean temp 

v2_out = nc.Dataset('/short/v45/lxy581/mom6/archive/so_mom6_v2/output009/prog.nc','r')
 
xh = v2_out.variables['xh'][:]
yh = v2_out.variables['yh'][:]

temp2 = v2_out.variables['temp'][:,0,:,:]
temp1 = np.nanmean(temp2,axis=0) 
temp = np.nanmean(temp1,axis=1)

plt.figure(2,figsize=(8,5))

plt.plot(temp,yh,'k')
plt.scatter(temp,yh,color='k')
plt.gca().xaxis.set_ticks_position('top')
plt.gca().xaxis.set_label_position('top')
plt.xlim(0,12)
plt.ylim(-1250,1250)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
#plt.xlabel('PTEMP ($^\circ$C)',fontsize=20)
plt.ylabel('Y (km)',fontsize=20)
plt.gca().set_position([0.2,0.1,0.6,0.8])


#plt.savefig('/short/v45/lxy581/mom6/diag/v2_temp.png',dpi=600)

plt.show()

