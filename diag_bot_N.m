% diag bottom stratification

disp('Running...');
clear;clc;

tic

addpath /home/581/lxy581/n2

dir0 = '/short/v45/lxy581/mom6/archive/so_mom6_v33/output000/';
dir  = '/short/v45/lxy581/mom6/archive/so_mom6_v33/output028/';

% read data
xh = ncread([dir,'prog.nc'],'xh');
yh = ncread([dir,'prog.nc'],'yh');
zl = ncread([dir,'prog.nc'],'zl');
t  = ncread([dir,'prog.nc'],'Time');

nx = length(xh);
ny = length(yh);
nz = length(zl);
nt = length(t);

% local f
f = ncread([dir0,'ocean_geometry.nc'],'f');

% local lat
DEG2RAD = pi/180;
OMEGA   = 7.292e-5;     %s-1   A.E.Gill p.597
% f     = 2*OMEGA*sin(lat*DEG2RAD);
lat     = asin(f/2/OMEGA)/DEG2RAD;

% reference density
rho = 1035;

% local depth
for j = 1:nx;
  for k = 1:ny;
    Z(j,k,:) = zl(:); % x,y,depth
  end
end 

for ii = 1:nz;
  LAT(:,:,ii) = lat(:,:);
end

% local g
g = sw_g(LAT,-Z);
n2m = zeros(nx,ny,nz-1); % x,y,z
num = 0;

for i = 1:nt;
  disp(i);

  % potential temperature
  temp = ncread([dir,'prog.nc'],'thetao',[1 1 1 i],[nx ny nz 1]);
  temp = squeeze(temp);

  % salinity
  salt = ncread([dir,'prog.nc'],'so',[1 1 1 i],[nx ny nz 1]);
  salt = squeeze(salt);

  % sea water pressure at sea water surface
  pso = ncread([dir,'prog.nc'],'pso',[1 1 i],[nx ny 1]);
  pso = squeeze(pso);

  % pressure
  for j = 1:nz;
    pres(:,:,j) = (pso(:,:) + rho*g(:,:,j)*zl(j))/1.e+4;
  end

  [l,n,m] = size(pres); % m-depth; n-y; l-x
  iup   = 1:m-1;
  ilo   = 2:m;
  p_ave = (pres(:,:,iup)+pres(:,:,ilo) )/2;
  pden_up = sw_pden(salt(:,:,iup),temp(:,:,iup),pres(:,:,iup),p_ave);
  pden_lo = sw_pden(salt(:,:,ilo),temp(:,:,ilo),pres(:,:,ilo),p_ave);

  mid_pden = (pden_up + pden_lo )/2;
  dif_pden = pden_up - pden_lo;
  mid_g    = (g(:,:,iup)+g(:,:,ilo))/2;
  dif_z    = diff(Z,1,3);
  n2       = -mid_g .* dif_pden ./ (dif_z .* mid_pden);      
  
  num = num + 1;
  n2m = n2m + n2;
end

n2m = n2m/num;

toc

tic

% find average values for the lowest 500m
disp('Finding and Averaging...')
b5up = ncread(['/short/v45/lxy581/mom6/diag/','find_b5.nc'],'b5up');
b5lo = ncread(['/short/v45/lxy581/mom6/diag/','find_b5.nc'],'b5low');

fac = zeros(nx,ny,nz-1);
for j = 1:ny;
  for i = 1:nx;
    fac(i,j,b5up(i,j)+1:b5lo(i,j)) = 1.0;
  end
end

n2m = n2m .* fac;
n2s = sum(n2m,3);

disp('Writing...')
%save /short/v45/lxy581/mom6/diag/n2_SO_10_b5_v32_y28_170825.mat n2m n2s xh yh
%save /short/v45/lxy581/mom6/diag/n2_SO_10_b5_v32_y30_170831.mat n2m n2s xh yh
save /short/v45/lxy581/mom6/diag/n2_SO_10_b5_v33_y50_170907.mat n2m n2s xh yh

toc
disp('Finished');
