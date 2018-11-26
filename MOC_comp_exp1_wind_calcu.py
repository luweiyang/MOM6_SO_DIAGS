import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cmocean

print 'Start...'

# Cd case, exp1, original topography, all wind cases
data_11 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1/output015/prog.nc','r')
data_12 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w2/output015/prog.nc','r')
data_13 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w3/output015/prog.nc','r')
data_11p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1p5/output015/prog.nc','r')
data_12p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w2p5/output015/prog.nc','r')
data_13p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w3p5/output015/prog.nc','r')

rhod_11 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1/output015/prog_rho.nc','r')
rhod_12 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w2/output015/prog_rho.nc','r')
rhod_13 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w3/output015/prog_rho.nc','r')
rhod_11p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w1p5/output015/prog_rho.nc','r')
rhod_12p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w2p5/output015/prog_rho.nc','r')
rhod_13p5 = nc.Dataset('/g/data1/v45/lxy581/payu/archive/exp1_w3p5/output015/prog_rho.nc','r')

coord = nc.Dataset('/short/v45/lxy581/mom6/input/exp1_w1/coord.nc','r')
rhol = coord.variables['Layer'][:]

yq = rhod_11.variables['yq'][:]
zl = rhod_11.variables['zl'][:]
t  = data_11.variables['Time'][:]
ny = np.size(yq)
nz = np.size(zl)
nt = np.size(t)
print 'nt = ',nt

# time-mean
vh_rho3_11 = np.nanmean(rhod_11.variables['vh_rho'][:,:,:,:],axis=0)
vh_rho3_12 = np.nanmean(rhod_12.variables['vh_rho'][:,:,:,:],axis=0)
vh_rho3_13 = np.nanmean(rhod_13.variables['vh_rho'][:,:,:,:],axis=0)
vh_rho3_11p5 = np.nanmean(rhod_11p5.variables['vh_rho'][:,:,:,:],axis=0)
vh_rho3_12p5 = np.nanmean(rhod_12p5.variables['vh_rho'][:,:,:,:],axis=0)
vh_rho3_13p5 = np.nanmean(rhod_13p5.variables['vh_rho'][:,:,:,:],axis=0)

# zonal-integral
vh_rho2_11 = np.nansum(vh_rho3_11,axis=-1)
vh_rho2_12 = np.nansum(vh_rho3_12,axis=-1)
vh_rho2_13 = np.nansum(vh_rho3_13,axis=-1)
vh_rho2_11p5 = np.nansum(vh_rho3_11p5,axis=-1)
vh_rho2_12p5 = np.nansum(vh_rho3_12p5,axis=-1)
vh_rho2_13p5 = np.nansum(vh_rho3_13p5,axis=-1)

# vertical cumulative sum from surface
vh_rhoc_11 = np.full((nz,ny),np.nan)
vh_rhoc_12 = np.full((nz,ny),np.nan)
vh_rhoc_13 = np.full((nz,ny),np.nan)
vh_rhoc_11p5 = np.full((nz,ny),np.nan)
vh_rhoc_12p5 = np.full((nz,ny),np.nan)
vh_rhoc_13p5 = np.full((nz,ny),np.nan)
vh_rhoc_11 = np.cumsum(vh_rho2_11,axis=0)
vh_rhoc_12 = np.cumsum(vh_rho2_12,axis=0)
vh_rhoc_13 = np.cumsum(vh_rho2_13,axis=0)
vh_rhoc_11p5 = np.cumsum(vh_rho2_11p5,axis=0)
vh_rhoc_12p5 = np.cumsum(vh_rho2_12p5,axis=0)
vh_rhoc_13p5 = np.cumsum(vh_rho2_13p5,axis=0)

# convert masked array to numpy array
vh_rhoc_11 = vh_rhoc_11.filled(np.nan)
vh_rhoc_12 = vh_rhoc_12.filled(np.nan)
vh_rhoc_13 = vh_rhoc_13.filled(np.nan)
vh_rhoc_11p5 = vh_rhoc_11p5.filled(np.nan)
vh_rhoc_12p5 = vh_rhoc_12p5.filled(np.nan)
vh_rhoc_13p5 = vh_rhoc_13p5.filled(np.nan)

# mask vh_rhoc for last few layers
for j in xrange(0,ny):
    for k in xrange(0,nz-1):
        if vh_rhoc_11[k+1,j] == vh_rhoc_11[k,j]:
            vh_rhoc_11[k+1,j] = 0.0
        if vh_rhoc_12[k+1,j] == vh_rhoc_12[k,j]:
            vh_rhoc_12[k+1,j] = 0.0
        if vh_rhoc_13[k+1,j] == vh_rhoc_13[k,j]:
            vh_rhoc_13[k+1,j] = 0.0
        if vh_rhoc_11p5[k+1,j] == vh_rhoc_11p5[k,j]:
            vh_rhoc_11p5[k+1,j] = 0.0
        if vh_rhoc_12p5[k+1,j] == vh_rhoc_12p5[k,j]:
            vh_rhoc_12p5[k+1,j] = 0.0
        if vh_rhoc_13p5[k+1,j] == vh_rhoc_13p5[k,j]:
            vh_rhoc_13p5[k+1,j] = 0.0

# mark the min & max boundary for surface density
rhont_11 = np.full((1,ny),np.nan)
rhomt_11 = np.full((1,ny),np.nan)
rhon_11 = np.zeros((1,ny))
rhom_11 = np.zeros((1,ny))

rhont_12 = np.full((1,ny),np.nan)
rhomt_12 = np.full((1,ny),np.nan)
rhon_12 = np.zeros((1,ny))
rhom_12 = np.zeros((1,ny))

rhont_13 = np.full((1,ny),np.nan)
rhomt_13 = np.full((1,ny),np.nan)
rhon_13 = np.zeros((1,ny))
rhom_13 = np.zeros((1,ny))

rhont_11p5 = np.full((1,ny),np.nan)
rhomt_11p5 = np.full((1,ny),np.nan)
rhon_11p5 = np.zeros((1,ny))
rhom_11p5 = np.zeros((1,ny))

rhont_12p5 = np.full((1,ny),np.nan)
rhomt_12p5 = np.full((1,ny),np.nan)
rhon_12p5 = np.zeros((1,ny))
rhom_12p5 = np.zeros((1,ny))

rhont_13p5 = np.full((1,ny),np.nan)
rhomt_13p5 = np.full((1,ny),np.nan)
rhon_13p5 = np.zeros((1,ny))
rhom_13p5 = np.zeros((1,ny))

for k in xrange(0,nt):
    print '-----'
    print k
    # for each time moment, min & max in zonal of surface density
    rho_11 = data_11.variables['rhoinsitu'][k,:,:,:]
    rhont_11[0,:] = np.nanmin(rho_11[0,:,:],axis=-1) 
    rhomt_11[0,:] = np.nanmax(rho_11[0,:,:],axis=-1)
    rhon_11 = np.append(rhon_11,rhont_11,0)
    rhom_11 = np.append(rhom_11,rhomt_11,0)
    rho_12 = data_12.variables['rhoinsitu'][k,:,:,:]
    rhont_12[0,:] = np.nanmin(rho_12[0,:,:],axis=-1) 
    rhomt_12[0,:] = np.nanmax(rho_12[0,:,:],axis=-1)
    rhon_12 = np.append(rhon_12,rhont_12,0)
    rhom_12 = np.append(rhom_12,rhomt_12,0)
    rho_13 = data_13.variables['rhoinsitu'][k,:,:,:]
    rhont_13[0,:] = np.nanmin(rho_13[0,:,:],axis=-1) 
    rhomt_13[0,:] = np.nanmax(rho_13[0,:,:],axis=-1)
    rhon_13 = np.append(rhon_13,rhont_13,0)
    rhom_13 = np.append(rhom_13,rhomt_13,0)
    rho_11p5 = data_11p5.variables['rhoinsitu'][k,:,:,:]
    rhont_11p5[0,:] = np.nanmin(rho_11p5[0,:,:],axis=-1) 
    rhomt_11p5[0,:] = np.nanmax(rho_11p5[0,:,:],axis=-1)
    rhon_11p5 = np.append(rhon_11p5,rhont_11p5,0)
    rhom_11p5 = np.append(rhom_11p5,rhomt_11p5,0)
    rho_12p5 = data_12p5.variables['rhoinsitu'][k,:,:,:]
    rhont_12p5[0,:] = np.nanmin(rho_12p5[0,:,:],axis=-1) 
    rhomt_12p5[0,:] = np.nanmax(rho_12p5[0,:,:],axis=-1)
    rhon_12p5 = np.append(rhon_12p5,rhont_12p5,0)
    rhom_12p5 = np.append(rhom_12p5,rhomt_12p5,0)
    rho_13p5 = data_13p5.variables['rhoinsitu'][k,:,:,:]
    rhont_13p5[0,:] = np.nanmin(rho_13p5[0,:,:],axis=-1) 
    rhomt_13p5[0,:] = np.nanmax(rho_13p5[0,:,:],axis=-1)
    rhon_13p5 = np.append(rhon_13p5,rhont_13p5,0)
    rhom_13p5 = np.append(rhom_13p5,rhomt_13p5,0)

print 'Loop finished'

rhon_11 = rhon_11[1:,:]
rhom_11 = rhom_11[1:,:]
rhon_11 = np.nanmin(rhon_11,axis=0)
rhom_11 = np.nanmax(rhom_11,axis=0)

rhon_12 = rhon_12[1:,:]
rhom_12 = rhom_12[1:,:]
rhon_12 = np.nanmin(rhon_12,axis=0)
rhom_12 = np.nanmax(rhom_12,axis=0)

rhon_13 = rhon_13[1:,:]
rhom_13 = rhom_13[1:,:]
rhon_13 = np.nanmin(rhon_13,axis=0)
rhom_13 = np.nanmax(rhom_13,axis=0)

rhon_11p5 = rhon_11p5[1:,:]
rhom_11p5 = rhom_11p5[1:,:]
rhon_11p5 = np.nanmin(rhon_11p5,axis=0)
rhom_11p5 = np.nanmax(rhom_11p5,axis=0)

rhon_12p5 = rhon_12p5[1:,:]
rhom_12p5 = rhom_12p5[1:,:]
rhon_12p5 = np.nanmin(rhon_12p5,axis=0)
rhom_12p5 = np.nanmax(rhom_12p5,axis=0)

rhon_13p5 = rhon_13p5[1:,:]
rhom_13p5 = rhom_13p5[1:,:]
rhon_13p5 = np.nanmin(rhon_13p5,axis=0)
rhom_13p5 = np.nanmax(rhom_13p5,axis=0)

print 'Saving data...'

# save data
fileobj = nc.Dataset('/short/v45/lxy581/mom6/diag/exp1_moc.nc','w')

fileobj.createDimension('y',ny)
fileobj.createDimension('level',nz)

yq_var = fileobj.createVariable('yq','f',('y',))
rho_var = fileobj.createVariable('rho','f',('level',))

var_11 = fileobj.createVariable('vh_11','f',('level','y'),fill_value=-1e+20)
var_12 = fileobj.createVariable('vh_12','f',('level','y'),fill_value=-1e+20)
var_13 = fileobj.createVariable('vh_13','f',('level','y'),fill_value=-1e+20)
var_11p5 = fileobj.createVariable('vh_11p5','f',('level','y'),fill_value=-1e+20)
var_12p5 = fileobj.createVariable('vh_12p5','f',('level','y'),fill_value=-1e+20)
var_13p5 = fileobj.createVariable('vh_13p5','f',('level','y'),fill_value=-1e+20)

var_min_11 = fileobj.createVariable('rhon_11','f',('y'),fill_value=-1e+20)
var_min_12 = fileobj.createVariable('rhon_12','f',('y'),fill_value=-1e+20)
var_min_13 = fileobj.createVariable('rhon_13','f',('y'),fill_value=-1e+20)
var_min_11p5 = fileobj.createVariable('rhon_11p5','f',('y'),fill_value=-1e+20)
var_min_12p5 = fileobj.createVariable('rhon_12p5','f',('y'),fill_value=-1e+20)
var_min_13p5 = fileobj.createVariable('rhon_13p5','f',('y'),fill_value=-1e+20)

var_max_11 = fileobj.createVariable('rhom_11','f',('y'),fill_value=-1e+20)
var_max_12 = fileobj.createVariable('rhom_12','f',('y'),fill_value=-1e+20)
var_max_13 = fileobj.createVariable('rhom_13','f',('y'),fill_value=-1e+20)
var_max_11p5 = fileobj.createVariable('rhom_11p5','f',('y'),fill_value=-1e+20)
var_max_12p5 = fileobj.createVariable('rhom_12p5','f',('y'),fill_value=-1e+20)
var_max_13p5 = fileobj.createVariable('rhom_13p5','f',('y'),fill_value=-1e+20)

yq_var[:] = yq[:]
rho_var[:] = rhol[:]

var_11[:,:] = vh_rhoc_11[:,:]
var_12[:,:] = vh_rhoc_12[:,:]
var_13[:,:] = vh_rhoc_13[:,:]
var_11p5[:,:] = vh_rhoc_11p5[:,:]
var_12p5[:,:] = vh_rhoc_12p5[:,:]
var_13p5[:,:] = vh_rhoc_13p5[:,:]

var_min_11[:] = rhon_11[:] 
var_min_12[:] = rhon_12[:] 
var_min_13[:] = rhon_13[:] 
var_min_11p5[:] = rhon_11p5[:] 
var_min_12p5[:] = rhon_12p5[:] 
var_min_13p5[:] = rhon_13p5[:] 

var_max_11[:] = rhom_11[:] 
var_max_12[:] = rhom_12[:] 
var_max_13[:] = rhom_13[:] 
var_max_11p5[:] = rhom_11p5[:] 
var_max_12p5[:] = rhom_12p5[:] 
var_max_13p5[:] = rhom_13p5[:] 

print 'Finsihed'
