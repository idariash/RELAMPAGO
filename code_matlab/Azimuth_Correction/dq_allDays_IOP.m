% Ivan Arias
% 2019/08/08

% azimuth_correction_mean = [];
% azimuth_correction_std = [];
% azimuth_correction_time = [];
% azimuth_correction_Allvalues =  [];
% azimuth_correction_AllvaluesTime = [];

%%



azimuth_correction_mean = [azimuth_correction_mean azimuthError_mean];
azimuth_correction_std = [azimuth_correction_std azimuthError_std];
azimuth_correction_time = [azimuth_correction_time azimuthError_time];
azimuth_correction_Allvalues =  [azimuth_correction_Allvalues azimuth_error_total];
azimuth_correction_AllvaluesTime = [azimuth_correction_AllvaluesTime scan_numTime];

%%
figure
%scatter(scan_time, azimuth_error_total)
hold on
scatter(azimuth_correction_AllvaluesTime, azimuth_correction_Allvalues, '.')
errorbar(azimuth_correction_time, azimuth_correction_mean, azimuth_correction_std)
datetick('x', 'HH:MM')
grid on 
%xlabel('HH:MM (UTC)')
datetick('x', 'mm/dd')
ylabel('Azimuth Difference (deg.)')
title(['Drift Evolution for IOP (2018)'])
hold off

%%
% a = azimuth_correction_mean';
% a = datetime(azimuth_correction_time,'ConvertFrom','datenum');
% a = a';
