% Ivan Arias
% Self Consistency
% CHIVO RELAMPAGO
% 2018/02/22

close all

filename = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';
filename_rvp9 = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/NetCDF/20181126/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';
filename_drops = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';


rhoHV_threshold = 0.98;

Z = ncread(filename, 'DBZ');
ZDR = ncread(filename, 'ZDR');
RHOHV = ncread(filename, 'RHOHV');
PHIDP = ncread(filename, 'PHIDP');
KDP = ncread(filename, 'KDP');
range = ncread(filename, 'range');

KDP_rvp9 = ncread(filename_rvp9, 'KDP');

Z(RHOHV < rhoHV_threshold) = NaN;
ZDR(RHOHV < rhoHV_threshold) = NaN;
VEL(RHOHV < rhoHV_threshold) = NaN;
PHIDP(RHOHV < rhoHV_threshold) = NaN;
KDP(RHOHV < rhoHV_threshold) = NaN;

KDP_indexing = (KDP_rvp9 > 0.009 & KDP_rvp9 < 0.021);
KDP(KDP_indexing) = nan;

%% Second sweep

Z = Z(360*length(range)+1:720*length(range));
ZDR = ZDR(360*length(range)+1:720*length(range));
RHOHV = RHOHV(360*length(range)+1:720*length(range));
PHIDP = PHIDP(360*length(range)+1:720*length(range));
KDP = KDP(360*length(range)+1:720*length(range));

Z(Z< 25) = nan;

%Estimation
Z_H = 10.^((Z)/10);
C = 1.46e-4;
alpha = 0.98;
beta = 0.2;
KDP_Estimated = C*(Z_H.^alpha).*10.^(-beta*ZDR)*180;
KDP_Estimated(KDP_Estimated<0.01) = nan;
KDP_Estimated(KDP_Estimated>5) = nan;
KDP(KDP<0) = nan;
close all
hist(KDP_Estimated,100)
figure
hist(KDP,100)
figure
Z_H_estimated =(1/C*KDP./(10.^(-beta*ZDR))).^(1/alpha);
Z_estimated = 10*log10(Z_H_estimated);

scatter(KDP, KDP_Estimated);
figure
scatter(Z, Z_Estimated);

return 
figure
%Z_estimated(Z_estimated<10) = nan;
%Z(Z<10) = nan;
scatter(Z, Z_estimated)
x = -10:1:60;
hold on 
plot(x,x)

figure
hist(Z_estimated,32)
figure
hist(Z,32)

Z_high = Z(Z > 40);
Z_high_estimated = real(Z_estimated(Z>40));
nanmean(Z_high - Z_high_estimated)

