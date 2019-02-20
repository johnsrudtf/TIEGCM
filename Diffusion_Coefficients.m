
%Author: Torfinn Johnsrud
%Created: Feb. 6 2019

clear all
% close all
clc

aa1=['/home/tojo5760/Documents/MATLAB/'];
linear=0;



%----------------
ut_want=1;% what time segment desired from simulation
feat = 2;%Select latitude and longitude desired
pdrag = 1;%1 if pdrag file used, 0 if not

if feat == 1
lon_want=100;%North nighttime feature
lat_want=63.75;%North nighttime feature
end
if feat == 2
lon_want=102.5;%South nighttime feature
lat_want=-58.75;%South nighttime feature
end
if feat == 3
lat_want=58.75;%Intermediate density area
lon_want=40;%Intermediate density area
end
if feat == 4 
lon_want=-150;% Daytime feature
lat_want=11.25;% Daytime feature
end
if feat == 5
lon_want=95;%Select Longitude desired
lat_want=62.5;%Select Latitude desired
end
%------------------
atom_unit=1.6605e-27; % kg/unit
%-----------------

%-----Loading Viki's tiegcm simulation-----
% format is (lon,lat,ilev,UT) format
if pdrag == 1
    filename = 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc';
    id = 'pdrag';
end
if pdrag == 0
    filename = 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc';
    id = 'no Ion Drag';
end
den = ncread(filename,'DEN');
zg = ncread(filename,'ZG'); %Geometric altitude
he = ncread(filename,'HE'); %Units of mass mixing ratio
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
% ar = ncread(filename,'AR');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32);%Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');


% fixed UT time 
den_tp1=squeeze(den(:,:,:,ut_want+1));% Total neutral density
zg_tp1=squeeze(zg(:,:,:,ut_want+1));% Geometric height
% z_tp1=squeeze(z(:,:,:,ut_want+1));% Geopotential height
he_tp1=squeeze(he(:,:,:,ut_want+1));% Helium mass mixing ratio
N2_tp1=squeeze(n2(:,:,:,ut_want+1));
O1_tp1=squeeze(o1(:,:,:,ut_want+1));
O2_tp1=squeeze(o2(:,:,:,ut_want+1));
% AR_tp1=squeeze(ar(:,:,:,ut_want+1));
tn_tp1=squeeze(tn(:,:,:,ut_want+1));% Neutral temperature
mbar_tp1=squeeze(mbar(:,:,:,ut_want+1));% Mean molecular mass
wn_tp1 = squeeze(wn(:,:,:,ut_want+1));


%fixed latitude and longitude  
i_index=find(lat==lat_want);

% Condense to selected latitude
den_alt_1=squeeze(den_tp1(:,i_index,:));
he_alt_1=squeeze(he_tp1(:,i_index,:));%TIEGCM output is mass mixing ratio
N2_alt_1=squeeze(N2_tp1(:,i_index,:));
O1_alt_1=squeeze(O1_tp1(:,i_index,:));
O2_alt_1=squeeze(O2_tp1(:,i_index,:));
% AR_alt_1=squeeze(Ar_tp1(:,i_index,:));
geometric_alt_1=squeeze(zg_tp1(:,i_index,:))/1e5; %Convert from cm to km
% geop_alt_1=squeeze(z_tp1(:,i_index,:))/1e5;
tn_alt_1=squeeze(tn_tp1(:,i_index,:));
mbar_alt_1=squeeze(mbar_tp1(:,i_index,:));
wn_alt_1 = squeeze(wn_tp1(:,i_index,:));

% -----Mean over the Lonitude-----
% den_alt=mean(den_alt_1);
% he_alt=mean(he_alt_1);
% N2_alt=mean(N2_alt_1);
% O1_alt=mean(O1_alt_1);
% % geop_alt=mean(geop_alt_1); 
% geom_alt=mean(geometric_alt_1);
% tn_alt=mean(tn_alt_1);
% mbar_alt=mean(mbar_alt_1);
% wn_alt = mean(wn_alt_1);

% -----Condense to selected Longitude-----
i_index = find(lon==lon_want);
den_alt = squeeze(den_alt_1(i_index,:));
he_alt = squeeze(he_alt_1(i_index,:));
N2_alt = squeeze(N2_alt_1(i_index,:));
O1_alt = squeeze(O1_alt_1(i_index,:));
O2_alt = squeeze(O2_alt_1(i_index,:));
% AR_alt = squeeze(AR_alt_1(i_index,:));
% geop_alt = squeeze(geop_alt_1(i_index,:));
geom_alt = squeeze(geometric_alt_1(i_index,:));
tn_alt = squeeze(tn_alt_1(i_index,:));
mbar_alt = squeeze(mbar_alt_1(i_index,:));%kg/kmol
wn_alt = squeeze(wn_alt_1(i_index,:));

% Number Density
nhe = he_alt.*den_alt*1e3/4/atom_unit;   %helium number density per m^3
nN2 = N2_alt.*den_alt*1e3/28.01/atom_unit; % N2 number density per m^3
nO1 = O1_alt.*den_alt*1e3/16/atom_unit; % O1 number density per m^3
nO2 = O2_alt.*den_alt*1e3/32/atom_unit; % O2 number density per m^3
% nAR = AR_alt.*den_alt*1e3/39.95/atom_unit; % Argon Number Density per m^3

%Mass Density
mhe = he_alt.*den_alt;% Helium Mass Density g/cm^3 
mN2 = N2_alt.*den_alt;% N2 Mass Density
mO1 = O1_alt.*den_alt;% Oxygen Mass Density


r=6372;%Earth Radius
g_tp=9.8.*r^2./(r+geom_alt).^2; %Find g as altitude changes
k=1.38064852e-23;% Boltzman's Constant
atom_unit=1.6605e-27;%kg/unit
Av = 6.022141*10^23;
kmol=6.02214*10^26;
mmw_he=0.004; %Helium atomic mass kg/mol
mmw_N2=0.02801; %N2 molecular mass
mmw_O1=0.016; % O1 molecular mass
mmw_O2=0.032; %O2 molecular mass
points = length(tn_alt);
meanmass = mbar_alt/1000;

%%-----Diffusion Calculation-----%%
Mass_O2 = 32.00;
Mass_He = 4.00;
Mass_N2 = 28.01;
Mass_Ar = 39.95;
Radius_O2 = 346/2.0;
Radius_He = 260/2.0;
Radius_N2 = 364/2.0;
Radius_Ar = 340/2.0;

Diff_denom_He = mmw_he*(collision_freq(nhe, nN2, Mass_He, Mass_N2, Radius_He, Radius_N2, tn_alt, tn_alt)+...
    collision_freq(nhe, nO2, Mass_He, Mass_O2, Radius_He, Radius_O2, tn_alt, tn_alt));
Diffusion_coeff_He = k*tn_alt./Diff_denom_He;

Diff_denom_O2 = mmw_O2*(collision_freq(nO2, nN2, Mass_O2, Mass_N2, Radius_O2, Radius_N2, tn_alt, tn_alt)+...
    collision_freq(nO2, nhe, Mass_O2, Mass_He, Radius_O2, Radius_He, tn_alt, tn_alt));
Diffusion_coeff_O2 = k*tn_alt./Diff_denom_O2;

Diff_denom_N2 = mmw_O2*(collision_freq(nN2, nO2, Mass_N2, Mass_O2, Radius_N2, Radius_O2, tn_alt, tn_alt)+...
    collision_freq(nN2, nhe, Mass_N2, Mass_He, Radius_N2, Radius_He, tn_alt, tn_alt));
Diffusion_coeff_N2 = k*tn_alt./Diff_denom_N2;

plot(Diffusion_coeff_He(1:end-12),geom_alt(1:end-12));
hold on
plot(Diffusion_coeff_O2(1:end-10),geom_alt(1:end-10));
plot(Diffusion_coeff_N2(1:end-10),geom_alt(1:end-10));
xlabel('Diffusion Coefficient');
ylabel('Altitude (km)');
title('Molecular Diffusion Coefficient');
legend('He','O2','N2');
