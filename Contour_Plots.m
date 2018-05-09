%Author: Torfinn Johnsrud
%Created: Fall 2017
%Updated: Feb 19. 2018
%Creates a latitude and longitude mesh to be contour plotted. Values can be any
% of the TIEGCM outputs. The mesh is saved as 2D matrix and written out to a
% text file. Normalized values or absolute values can be output.
% I used Python to create contour plots because I liked the way they looked more.
% Note, files are assumed to be finer TIEGCM output with 2.5 degrees
% between data points.

clear all
close all
clc

aa1=['/home/tojo5760/Documents/TIEGCM_files/'];
linear=0;



%----------------
ut_want=1; % Specify the universal time. This is the first UT time in the TIEGCM file
alt_want=400; %Specify the desired altitude in km
pdrag = 0;% set to 1 to use pdrag file, 0 to use ctrSS file

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

%Read in data arrays
den = ncread(filename,'DEN');%g/cm^3
zg = ncread(filename,'ZG')/1e5;%geometric height in km
he = ncread(filename,'HE'); %Units of TIEGCM's mass mixing ratio (eg. mass of species/total mass)
n2 = ncread(filename,'N2');
o1 = ncread(filename,'O1');
o2 = ncread(filename,'O2');
tn = ncread(filename,'TN'); %Neutral Temperature
mbar = 1./(he/4+n2/28+o1/16+o2/32); %Getting mean molecular mass from mmr
lat = ncread(filename,'lat');
lon = ncread(filename,'lon');
wn = ncread(filename,'WN'); %Neutral Winds


% fixed UT time. Condense arrays to desired UT time. Note ut_want+1 is used due to errors with ut_want
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
%Performs a linear interpolation to get values at the desired altitude

geom_index_1 = zg-alt_want;

for n=1:72 %Number of sample points for latitude (180/2.5)
    for m=1:144
        geom_index=geom_index_1(m,n,:);
        geom_index=squeeze(geom_index);
        [val, i_index] = min(abs(geom_index));%Find altitude closest to desired altitude
        if zg(m,n,i_index)>=alt_want %Linearly interpolate backwards
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
        if zg(m,n,i_index)<alt_want %Linearly interpolate forwards
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
left = Tn(:,(73:end)); %TIEGCM output is in longitudes of -180 to 180. This changes it to be 0 to 360
right = Tn(:,(1:72));
Tn = [left,right];
maxTn = max(max(Tn)); %Calculated for normalization
minTn = min(min(Tn));
rangeTn = maxTn-minTn;
mid_Tn = (maxTn+minTn)/2;
Tn_normalized = (Tn-mid_Tn)./rangeTn*2;%Scale data to between -1 and 1

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
mHelog = log10(mHe);% Misc. calculated values

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

%Shift N2 for max correlation----
%This finds the shift of N2 plot to be maximally correlated with helium
lowest = 1;
for i=2:72
    bottom = N2_normalized((i:end),:);
    top = N2_normalized((1:i-1),:);
    shiftvert_N2 = [top;bottom];% This shifts the N2 data up by (72-i)
    for j=2:144
        left = shiftvert_N2(:,(j:end));
        right = shiftvert_N2(:,(1:j-1));
        normN2_shift = [left,right];% Shifts N2 data to the right by (144-j)
        R = corrcoef(He_normalized, normN2_shift);
        correlation = R(1,2);
        if correlation<lowest
            lowest = correlation;
            location = [i,j];%
        end
    end
end

%---------------------------------
mO1 = O1.*Den*1000;% Density kg/m^3
mO1 = mO1.';
left = mO1(:,(73:end));
right = mO1(:,(1:72));
mO1 = [left,right];

ratio = mHe./mN2;
oxy_ratio = mO1./mN2;
normal_he_n2 = He_normalized./N2_normalized;
normal_he_tn = He_normalized./Tn_normalized;

Difference = He_normalized-Tn_normalized;

dlmwrite(['He_',id,'400km.txt'],mHe);

%Correlation coefficients
CorrN2 = corrcoef(He_normalized, N2_normalized);
CorrTn = normxcorr2(He_normalized, Tn_normalized);
Corr = normxcorr2(He_normalized, He_normalized);
%dlmwrite('He_N2_Correlation_400km_pdrag',CorrN2); %Write output to text
%file for plotting
