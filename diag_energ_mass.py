
# Diagnose and plot En, APE and KE from the beginning

# v2   output000 - output009 - output014

# v24  output000 - output009 - output014

# v32  output000 - output023

# Energy/Mass, Mass time series

import numpy as np 
import netCDF4 as nc
import matplotlib.pyplot as plt

En_t = np.zeros((1))
M_t  = np.zeros((1))
time = np.zeros((1))

for opt in xrange(0,24,1):
  print '--------'
  print opt
  print '/short/v45/lxy581/mom6/archive/so_mom6_v32/output0%02d/ocean.stats.nc' % opt
  dir_opt = '/short/v45/lxy581/mom6/archive/so_mom6_v32/output0%02d/ocean.stats.nc' % opt

  data = nc.Dataset(dir_opt,'r')

  t = data.variables['Time'][:]
  print np.shape(t)
 
  En = data.variables['En'][:]
  M  = data.variables['Mass'][:]

  time = np.append(time,t,0)
  En_t = np.append(En_t,En,0)
  M_t  = np.append(M_t,M,0)

time = time[1:]
En_t = En_t[1:]
M_t  = M_t[1:]
E_M  = En_t/M_t

# Plot
fig1 = plt.figure(1,figsize=(12,4))

ax1 = fig1.add_subplot(111)
l1 = ax1.plot(time,E_M,'k',linewidth=2.0,label='Energy/Mass')
ax1.grid(True)
ax1.set_xticks(np.arange(0,365*30+365*5,365*5))
ax1.set_xticklabels(['0','5','10','15','20','25','30'])
#ax1.set_xticks(np.arange(0,365*50+365*10,365*10))
#ax1.set_xticklabels(['0','10','20','30','40','50'])
ax1.tick_params(axis='x',labelsize=16)
ax1.tick_params(axis='y',labelsize=16)
ax1.set_xlabel('Year',fontsize=20)
ax1.set_yticks(np.arange(0,0.8e-2 + 0.2e-2,0.2e-2))
ax1.set_yticklabels(['0','0.2','0.4','0.6','0.8'])
#ax1.set_yticks(np.arange(0,4.0e-2,0.8e-2))
#ax1.set_yticklabels(['0','0.8','1.6','2.4','3.2'])
ax1.set_ylabel('Energy/Mass ($10^{-2}$ J)',fontsize=20)
ax1.set_position([0.1,0.2,0.8,0.6])

ax2 = fig1.add_subplot(111,sharex=ax1,frameon=False)
l3 = ax2.plot(time,M_t,'b',linewidth=2.0,label='Total Mass')
ax2.yaxis.tick_right()
ax2.tick_params(axis='x',labelsize=16)
ax2.tick_params(axis='y',labelsize=16,colors='blue')
ax2.yaxis.set_label_position("right")
#ax2.set_yticks(np.arange(0,0.12e+18,0.02e+18))
#ax2.set_yticklabels(['0','0.02','0.04','0.06','0.08','0.10'])
#ax2.set_yticks(np.arange(0,1.8e+17,0.3e+17))
#ax2.set_yticklabels(['0','0.3','0.6','0.9','1.2','1.5'])
ax2.set_yticks(np.arange(0,5.0e+19,1.0e+19))
ax2.set_yticklabels(['0','1.0','2.0','3.0','4.0'])
ax2.set_ylabel('Mass ($10^{19}$ J)',fontsize=20)
ax2.yaxis.label.set_color('blue')
ax2.set_position([0.1,0.2,0.8,0.6])

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=[0.78,0])

plt.savefig('/short/v45/lxy581/mom6/diag/v32_energ_mass.png',dpi=600)

plt.show()

