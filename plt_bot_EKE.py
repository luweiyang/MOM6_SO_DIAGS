# plot bot 500m EKE
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import scipy.io
#from mpl_toolkits.basemap import basemap

data = nc.Dataset('/short/v45/lxy581/mom6/diag/bot_500_SO_10_EKE_spd_170823.nc','r')

xh = data.variables['xh'][:]
yh = data.variables['yh'][:]
EKE = data.variables['mEKEb5'][:]
spd = data.variables['mspdb5'][:]

# plot - v32, output023, bot 500m-averaged speed

plt.figure(1,figsize=(8,5))

s_level = np.arange(0,0.2+0.01,0.01)
s_ticks = np.arange(0,0.2+0.05,0.05)

X,Y=np.meshgrid(xh,yh)

pc = plt.contourf(X,Y,spd,cmap=plt.cm.jet,levels=s_level,extend='both')
cb = plt.colorbar(pc,ticks=s_ticks)
plt.gca().set_yticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_yticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_xticks(np.arange(-2000,2000 + 1000,1000))
plt.gca().set_xticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.15,0.2,0.65,0.6])
cb.ax.set_position([0.85,0.2,0.03,0.6])
cb.set_label('Bottom Speed (m/s)',y=0.5,fontsize=16)
cb.ax.tick_params(labelsize=16)
plt.xlabel('X (km)',fontsize=16)
plt.ylabel('Y (km)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.savefig('/short/v45/lxy581/mom6/diag/v32_y28_b5_spd_170823.png',dpi=600)

plt.show()

