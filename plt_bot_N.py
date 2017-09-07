import numpy as np
import matplotlib.pyplot as plt
import scipy.io 

mat = scipy.io.loadmat('/short/v45/lxy581/mom6/diag/n2_SO_10_b5_v33_y50_170907.mat')

xh  = mat['xh']
yh  = mat['yh']
n2s = mat['n2s']
ns0 = np.sqrt(n2s)
ns0 = np.real(ns0)
ns  = np.log10(ns0)
ns  = np.transpose(ns)

plt.figure(1,figsize=(8,5))

n_level = np.arange(-4,-2+0.01,0.001) 
n_ticks = np.arange(-4,-2+1,1) 

X,Y=np.meshgrid(xh,yh)

pc = plt.contourf(X,Y,ns,cmap=plt.cm.jet,levels=n_level,extend='both')
cb = plt.colorbar(pc,ticks=n_ticks)
plt.gca().set_yticks(np.arange(-1250,1250 + 500,500))
plt.gca().set_yticklabels(['0','500','1000','1500','2000','2500'])
plt.gca().set_xticks(np.arange(-2000,2000 + 1000,1000))
plt.gca().set_xticklabels(['0','1000','2000','3000','4000'])
plt.gca().set_position([0.15,0.2,0.65,0.6])
cb.ax.set_position([0.85,0.2,0.03,0.6])
cb.set_label('Bottom Stratification ($log_{10}(s^{-1})$)',y=0.5,fontsize=16)
cb.ax.tick_params(labelsize=16)
plt.xlabel('X (km)',fontsize=16)
plt.ylabel('Y (km)',fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

#plt.savefig('/short/v45/lxy581/mom6/diag/v32_y28_b5_N_170825.png',dpi=600)
#plt.savefig('/short/v45/lxy581/mom6/diag/v32_y30_b5_N_170831.png',dpi=600)
plt.savefig('/short/v45/lxy581/mom6/diag/v33_y50_b5_N_170907.png',dpi=600)

plt.show()
