%Author: Torfinn Johnsrud
%Created: Fall 2017
%Updated: Feb 19. 2018
%Creates a latitude and longitude mesh to be contour plotted.

clear all
close all
clc

aa1=['/home/tojo5760/Documents/TIEGCM_files/'];
linear=1;



%----------------
ut_want=1;%
alt_want=400;
pdrag = 1;

%------------------
atom_unit=1.67e-27; % kg/unit
%-----------------

%-----Loading Viki's tiegcm simulation-----
% Follows (lon,lat,ilev,UT) format
if pdrag == 1
    filename = 'HSUVW.tiegcm2.0_dres.pdrag_f107_180_001.nc';
    id = 'pdrag';
end
if pdrag == 0
    filename = 'HSUVW.tiegcm2.0_dres.nodragtest_ctrSS_f107_180_001.nc';
    id = 'ctrSS';
end
den = ncread(filename,'DEN');%g/cm^3
zg = ncread(filename,'ZG')/1e5;%geometric height in km
he = ncread(filename,'HE');
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN');
mbar = 1./(he/4+n2/28+o1/16+o2/32);%Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN');


% fixed UT time 
den=squeeze(den(:,:,:,ut_want+1));% Total neutral density
zg=squeeze(zg(:,:,:,ut_want+1));% Geometric height
% z_tp1=squeeze(z(:,:,:,ut_want+1));% Geopotential height
he=squeeze(he(:,:,:,ut_want+1));% Helium mass mixing ratio
n2=squeeze(n2(:,:,:,ut_want+1));
o1=squeeze(o1(:,:,:,ut_want+1));
tn=squeeze(tn(:,:,:,ut_want+1));% Neutral temperature
mbar=squeeze(mbar(:,:,:,ut_want+1));% Mean molecular mass
wn = squeeze(wn(:,:,:,ut_want+1));


% -----Select Altitude-----

geom_index_1 = zg-alt_want;

for n=1:72 %Number of sample points for latitude (180/2.5)
    for m=1:144
        geom_index=geom_index_1(m,n,:);
        geom_index=squeeze(geom_index);
        [val, i_index] = min(abs(geom_index));%Find altitude closest to desired altitude
        if zg(m,n,i_index)>=alt_want %Linearly interpolate
            y0 = den(m,n,i_index-1);
            y1 = den(m,n,i_index);
            x0 = zg(m,n,i_index-1);
            x1 = zg(m,n,i_index);
            Den(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);

            y0 = tn(m,n,i_index-1);
            y1 = tn(m,n,i_index);
            Tn(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);

            y0 = he(m,n,i_index-1);
            y1 = he(m,n,i_index);
            He(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
            
            y0 = n2(m,n,i_index-1);
            y1 = n2(m,n,i_index);
            N2(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
            
            y0 = o1(m,n,i_index-1);
            y1 = o1(m,n,i_index);
            O1(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
        end
        if zg(m,n,i_index)<alt_want %Linearly interpolate
            y0 = den(m,n,i_index);
            y1 = den(m,n,i_index+1);
            x0 = zg(m,n,i_index);
            x1 = zg(m,n,i_index+1);
            Den(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);

            y0 = tn(m,n,i_index);
            y1 = tn(m,n,i_index+1);
            Tn(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);

            y0 = he(m,n,i_index);
            y1 = he(m,n,i_index+1);
            He(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
            
            y0 = n2(m,n,i_index);
            y1 = n2(m,n,i_index+1);
            N2(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
            
            y0 = o1(m,n,i_index);
            y1 = o1(m,n,i_index+1);
            O1(m,n) = (y0*(x1-alt_want)+y1*(alt_want-x0))/(x1-x0);
    end
end
end

%Temperature
Tn = Tn.';
left = Tn(:,(73:end));
right = Tn(:,(1:72));
Tn = [left,right];
maxTn = max(max(Tn));
minTn = min(min(Tn));
rangeTn = maxTn-minTn;
Tn_normalized = (Tn-minTn)./rangeTn;

%Mass Density
mHe = He.*Den*1000;% Helium Mass Density kg/m^3 
mHe = mHe.';
left = mHe(:,(73:end));
right = mHe(:,(1:72));
mHe = [left,right];
maxHe = max(max(mHe));
minHe = min(min(mHe));
rangeHe = maxHe-minHe;
mid_He = (maxHe+minHe)/2;
He_normalized = (mHe-mid_He)./rangeHe*2;%Scale data to between -1 and 1
mHelog = log10(mHe);

mN2 = N2.*Den*1000;%kg/m^3
mN2 = mN2.';
left = mN2(:,(73:end));
right = mN2(:,(1:72));
mN2 = [left,right];
maxN2 = max(max(mN2));
minN2 = min(min(mN2));
rangeN2 = maxN2-minN2;
mid_N2 = (maxN2+minN2)/2;
N2_normalized = (mN2-mid_N2)./rangeN2*2;
mN2log = log10(mN2);

mO1 = O1.*Den*1000;% Density kg/m^3 
mO1 = mO1.';
left = mO1(:,(73:end));
right = mO1(:,(1:72));
mO1 = [left,right];

ratio = mHe./mN2;
oxy_ratio = mO1./mN2;
normal_he_n2 = He_normalized./N2_normalized;
normal_he_tn = He_normalized./Tn_normalized;

dlmwrite(['He_400km_',id,'_TEST_DELETE.txt'],He_normalized);

%Correlation
CorrN2 = xcorr2(He_normalized, N2_normalized);
CorrTn = xcorr2(He_normalized, Tn_normalized);
Corr = xcorr2(He_normalized);
%dlmwrite('He_N2_Correlation_400km_pdrag',CorrN2);