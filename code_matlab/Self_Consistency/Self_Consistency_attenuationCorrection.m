% Ivan Arias
% Self Consistency
% CHIVO RELAMPAGO after attenuation correction using diffent coeficients
% 2018/10/02

addpath /net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab

filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/corcsapr2cfrppiM1.a1.20190126.053003_Nesbitt.nc';

%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/corcsapr2cfrppiM1.a1.20190126.053003_Book.nc';

%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR_Book.nc';

%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR_Nesbitt.nc';
%'/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';
filename_rvp9 = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/NetCDF/20181126/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';
filename_drops = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/cfrad.20181126_071002.505_to_20181126_071536.425_col-radar_PPICOWH360_SUR.nc';


rhoHV_threshold = 0.95;

Z = ncread(filename, 'corrected_reflectivity');
ZDR = ncread(filename, 'corrected_differential_reflectivity') - (-3.6);
RHOHV = ncread(filename, 'corrected_cross_correlation_ratio');
PHIDP = ncread(filename, 'corrected_differential_phase');
KDP = ncread(filename, 'corrected_specific_differential_phase');
range = ncread(filename, 'range');

%KDP_rvp9 = ncread(filename_rvp9, 'KDP');

Z(RHOHV < rhoHV_threshold) = NaN;
ZDR(RHOHV < rhoHV_threshold) = NaN;
VEL(RHOHV < rhoHV_threshold) = NaN;
PHIDP(RHOHV < rhoHV_threshold) = NaN;
KDP(RHOHV < rhoHV_threshold) = NaN;

%KDP_indexing = (KDP_rvp9 > 0.009 & KDP_rvp9 < 0.021);
%KDP(KDP_indexing) = nan;

%% Second sweep

Z = Z(360*length(range)+1:720*length(range));
ZDR = ZDR(360*length(range)+1:720*length(range));
RHOHV = RHOHV(360*length(range)+1:720*length(range));
PHIDP = PHIDP(360*length(range)+1:720*length(range));
KDP = KDP(360*length(range)+1:720*length(range));

Z(Z< 35) = nan;

%Estimation
Z_H = 10.^((Z)/10);
C = 1.46e-4;
alpha = 0.98;
beta = 0.2;
KDP_Estimated = C*(Z_H.^alpha).*10.^(-beta*ZDR)*180;
KDP_Estimated(KDP_Estimated<0.01) = nan;
KDP_Estimated(KDP_Estimated>5) = nan;
KDP(KDP<0) = nan;
%close all
% hist(KDP_Estimated,100)
% figure
% hist(KDP,100)
% figure
Z_H_estimated =(1/C*KDP./(10.^(-beta*ZDR))).^(1/alpha);
Z_estimated = 10*log10(Z_H_estimated);

%%
figure
[A,CZplot,CZDRplot]=p_icp(Z, Z_estimated,50);
mean_A = nanmean(A,2);
figure; pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
caxis([0 3]);
hold on 
plot(x,x, 'r')
hc=colorbar;
xlabel(hc,'log_{10}(N_{obs})');
xlim([25 60])
ylim([10 60])

Z_estimated_strong = Z_estimated;
Z_estimated_strong(Z_estimated_strong < 35) = nan;

rho = nanmean((Z - nanmean(Z)).*(Z_estimated_strong - nanmean(Z_estimated_strong)))/(nanstd(Z)*nanstd(Z_estimated_strong));
Bias = nanmean(Z - Z_estimated_strong);
RMSE = sqrt(nanmean((Z - Z_estimated_strong).^2));

title(['Bias: ' num2str(round(Bias,2)) 'dB | CORR: ' num2str(round(rho,2))...
    ' | RMSE: ' num2str(round(RMSE,2)) 'dB']);

return 
xlim([27 37])
ylim([27 37])
title(['Bias: ' num2str(round(Bias,2)) 'dB | CORR: ' num2str(round(rho,2))...
    ' | RMSE: ' num2str(round(RMSE,2)) 'dB']);


figure
scatter(KDP, KDP_Estimated);
figure
scatter(Z, Z_estimated);
hold off

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

