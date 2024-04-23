% Ivan Arias
% 2020/01/14
% Compute histrograms of the region between CSAPR and CHIVO

addpath /net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Codes/utiles

filename_chivo = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/DualRadar/20181130/DROPS/cfrad.20181130_033054.697_to_20181130_033709.927_col-radar_REL_PNL360A_SUR.nc';
filename_csapr = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20181130/DROPS/corcsapr2cfrppiM1.a1.20181130.033003.nc';
azimuth_offset = 3.06; % for 2018/11/30

DBZ_chivo = ncread(filename_chivo, 'corrected_reflectivity');
ZDR_chivo = ncread(filename_chivo, 'corrected_differential_reflectivity');
RHOHV_chivo =  ncread(filename_chivo, 'corrected_cross_correlation_ratio');
KDP_chivo = ncread(filename_chivo, 'corrected_specific_differential_phase');
range_chivo = double(ncread(filename_chivo, 'range'))'/1e3;
azimuth_chivo = ncread(filename_chivo, 'azimuth') - azimuth_offset;
elevation_chivo = ncread(filename_chivo, 'elevation');
ray_n_gates_chivo = ncread(filename_chivo, 'ray_n_gates');

[Range_chivo, Azimuth_chivo, Elevation_chivo] = get_range_azimuth_index(range_chivo, azimuth_chivo, elevation_chivo, ray_n_gates_chivo);
Height_chivo = Range_chivo.*sind(Elevation_chivo) + 0.460; % Height MSL
r_chivo = Range_chivo.*cosd(Elevation_chivo);
X_chivo = r_chivo.*cosd(90 - Azimuth_chivo);
Y_chivo = r_chivo.*sind(90 - Azimuth_chivo);

DBZ_csapr = ncread(filename_csapr, 'corrected_reflectivity');
ZDR_csapr = ncread(filename_csapr, 'corrected_differential_reflectivity');
RHOHV_csapr =  ncread(filename_csapr, 'corrected_cross_correlation_ratio');
KDP_csapr = ncread(filename_csapr, 'corrected_specific_differential_phase');
range_csapr = double(ncread(filename_csapr, 'range'))'/1e3;
azimuth_csapr = ncread(filename_csapr, 'azimuth') - azimuth_offset;
elevation_csapr = ncread(filename_csapr, 'elevation');
ray_n_gates_csapr = ncread(filename_csapr, 'ray_n_gates');

[Range_csapr, Azimuth_csapr, Elevation_csapr] = get_range_azimuth_index(range_csapr, azimuth_csapr, elevation_csapr, ray_n_gates_csapr);
Height_csapr = Range_csapr.*sind(Elevation_csapr) + 0.996; % Height MSL
r_csapr = Range_csapr.*cosd(Elevation_csapr);
X_csapr = r_csapr.*cosd(90 - Azimuth_csapr);
Y_csapr = r_csapr.*sind(90 - Azimuth_csapr);

% Isolating region beween the two radars

Index_chivo = abs(X_chivo - -54/2) < 8 & abs(Y_chivo - -54/2) < 8 & ... % CHIVO is in the III quadrant of the cartisian plane 
    180 < Azimuth_chivo & Azimuth_chivo < 270 & 1.7 < Height_chivo & Height_chivo < 2.8;
Index_csapr = abs(X_csapr - 54/2) < 8 & abs(Y_csapr - 54/2) < 8 & ...
    0 < Azimuth_csapr & Azimuth_csapr < 90 & 1.7 < Height_csapr & Height_csapr < 2.8;

DBZ_chivo(RHOHV_chivo < 0.95) = nan;
ZDR_chivo(RHOHV_chivo < 0.95) = nan;
dbz_chivo = DBZ_chivo(Index_chivo);
zdr_chivo = ZDR_chivo(Index_chivo);

DBZ_csapr(RHOHV_csapr < 0.95) = nan;
ZDR_csapr(RHOHV_csapr < 0.95) = nan;
dbz_csapr = DBZ_csapr(Index_csapr);
zdr_csapr = ZDR_csapr(Index_csapr);

%%
mean_chivo = nanmean(dbz_chivo);
mean_csapr = nanmean(dbz_csapr);
median_chivo = nanmedian(dbz_chivo);
median_csapr = nanmedian(dbz_csapr);
std_chivo = nanstd(dbz_chivo);
std_csapr = nanstd(dbz_csapr);
figure
hist(dbz_chivo, 32)
xlim([10, 55])
xlabel('Reflectivity (dBZ)')
ylabel('Counts')
title(['mean: ' num2str(round(mean_chivo, 2)) ' | median: ' num2str(round(median_chivo, 2))...
    ' | std: ' num2str(round(std_chivo, 2))])
grid on

figure
hist(dbz_csapr, 32)
xlim([10, 55])
xlabel('Reflectivity (dBZ)')
ylabel('Counts')
title(['mean: ' num2str(round(mean_csapr, 2)) ' | median: ' num2str(round(median_csapr, 2))...
    ' | std: ' num2str(round(std_csapr, 2))])
grid on

%%
mean_chivo = nanmean(zdr_chivo);
mean_csapr = nanmean(zdr_csapr);
median_chivo = nanmedian(zdr_chivo);
median_csapr = nanmedian(zdr_csapr);
std_chivo = nanstd(zdr_chivo);
std_csapr = nanstd(zdr_csapr);
figure
hist(zdr_chivo, 32)
xlim([-2, 4])
xlabel('Differential Reflectivity (dB)')
ylabel('Counts')
title(['mean: ' num2str(round(mean_chivo, 2)) ' | median: ' num2str(round(median_chivo, 2))...
    ' | std: ' num2str(round(std_chivo, 2))])
grid on

figure
hist(zdr_csapr, 32)
xlim([-2, 4])
xlabel('Differential Reflectivity (dB)')
ylabel('Counts')
title(['mean: ' num2str(round(mean_csapr, 2)) ' | median: ' num2str(round(median_csapr, 2))...
    ' | std: ' num2str(round(std_csapr, 2))])
grid on



