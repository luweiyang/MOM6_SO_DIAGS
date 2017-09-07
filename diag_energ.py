# Diagnose and plot En, APE and KE from the beginning

# v2   output000 - output009 - output014

# v32  output000 - output023

# v33  output000 - output028

# READ En, APE and KE from outputs


import numpy as np 
import netCDF4 as nc
import matplotlib.pyplot as plt

En_t  = np.zeros((1))
APE_t = np.zeros((1))
KE_t  = np.zeros((1))
M_t  = np.zeros((1))
time  = np.zeros((1))

for opt in xrange(0,28+1,1):
  print '--------'
  print opt
  print '/short/v45/lxy581/mom6/archive/so_mom6_v33/output0%02d/ocean.stats.nc' % opt
  dir_opt = '/short/v45/lxy581/mom6/archive/so_mom6_v33/output0%02d/ocean.stats.nc' % opt

  data = nc.Dataset(dir_opt,'r')

  t = data.variables['Time'][:]
  print np.shape(t)
 
  En  = data.variables['En'][:]
  APE = data.variables['APE'][:,:]
  KE  = data.variables['KE'][:,:]
  M   = data.variables['Mass'][:]

  # Total
  APE_sum = np.nansum(APE,axis=1)
  KE_sum  = np.nansum(KE,axis=1)
  
  time  = np.append(time,t,0)
  En_t  = np.append(En_t,En,0)
  M_t   = np.append(M_t,M,0)
  APE_t = np.append(APE_t,APE_sum,0)
  KE_t  = np.append(KE_t,KE_sum,0)


time  = time[1:]
En_t  = En_t[1:]
M_t   = M_t[1:]
APE_t = APE_t[1:]
KE_t  = KE_t[1:]

En_M  = En_t/M_t
APE_M = APE_t/M_t
KE_M  = KE_t/M_t

# Plot
fig1 = plt.figure(1,figsize=(12,4))

ax1 = fig1.add_subplot(111)
l1 = ax1.plot(time,En_M,'k',linewidth=2.0,label='Total Energy/Mass')
l2 = ax1.plot(time,KE_M,'k--',linewidth=1.0,label='Total KE/Mass')
ax1.grid(True)
ax1.set_xticks(np.arange(0,365*50+365*10,365*10))
ax1.set_xticklabels(['0','10','20','30','40','50'])
ax1.tick_params(axis='x',labelsize=16)
ax1.tick_params(axis='y',labelsize=16)
ax1.set_xlabel('Year',fontsize=20)
ax1.set_yticks(np.arange(0,4e-3 + 1e-3,1e-3))
ax1.set_yticklabels(['0','1','2','3','4'])
ax1.set_ylabel('Energy/Mass ($10^{-3}$ J/kg)',fontsize=20)
ax1.set_position([0.1,0.2,0.8,0.6])

ax2 = fig1.add_subplot(111,sharex=ax1,frameon=False)
l3 = ax2.plot(time,APE_M,'b',linewidth=2.0,label='Total APE/Mass')
ax2.yaxis.tick_right()
ax2.tick_params(axis='x',labelsize=16)
ax2.tick_params(axis='y',labelsize=16,colors='blue')
ax2.yaxis.set_label_position("right")
ax2.set_yticks(np.arange(0,6e-4 + 1.5e-4,1.5e-4))
ax2.set_yticklabels(['0','1.5','3.0','4.5','6.0'])
ax2.set_ylabel('Energy/Mass ($10^{-4}$ J/kg)',fontsize=20)
ax2.yaxis.label.set_color('blue')
ax2.set_position([0.1,0.2,0.8,0.6])

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=[0.68,0])

plt.savefig('/short/v45/lxy581/mom6/diag/v33_energ_mass.png',dpi=600)

plt.show()

