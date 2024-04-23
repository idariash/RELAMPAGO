% Ivan Arias
% 2019/09/30
load('IOP_drift.mat')
all_azimuth_correction_mean = azimuth_correction_mean;
all_azimuth_correction_std = azimuth_correction_std;
all_azimuth_correction_time = azimuth_correction_time;

load('ExtendedPeriod_drift.mat')
all_azimuth_correction_mean = [all_azimuth_correction_mean azimuth_correction_mean];
all_azimuth_correction_std = [all_azimuth_correction_std azimuth_correction_std];
all_azimuth_correction_time = [all_azimuth_correction_time azimuth_correction_time];

all_azimuth_correction_mean = all_azimuth_correction_mean';
all_azimuth_correction_std = all_azimuth_correction_std';
all_azimuth_correction_time = all_azimuth_correction_time';

all_azimuth_correction_mean = all_azimuth_correction_mean(~isnan(all_azimuth_correction_mean));
all_azimuth_correction_std = all_azimuth_correction_std(~isnan(all_azimuth_correction_mean));
all_azimuth_correction_time = all_azimuth_correction_time(~isnan(all_azimuth_correction_mean));

%figure
plot(all_azimuth_correction_time, all_azimuth_correction_mean)
datetick('x', 'HH:MM')
grid on 
%xlabel('HH:MM (UTC)')
datetick('x', 'mm/dd')
ylabel('Azimuth Difference (deg.)')
title('Drift Evolution for RELAMPAGO 2018 - 2019')

%%
all_mean_azimuth_correction_time = reshape(all_azimuth_correction_time, 4, 55);
all_mean_azimuth_correction_std = reshape(all_azimuth_correction_std, 4, 55);
all_mean_azimuth_correction_mean = reshape(all_azimuth_correction_mean, 4, 55);

all_mean_azimuth_correction_time = nanmean(all_mean_azimuth_correction_time);
all_mean_azimuth_correction_std = nanmean(all_mean_azimuth_correction_std)/2;
all_mean_azimuth_correction_mean = nanmean(all_mean_azimuth_correction_mean);


errorbar(all_mean_azimuth_correction_time, all_mean_azimuth_correction_mean, all_mean_azimuth_correction_std)
datetick('x', 'HH:MM')
grid on 
%xlabel('HH:MM (UTC)')
datetick('x', 'mm/dd')
ylabel('Mean and Std. of Azimuth Difference (deg.)')
title(['Drift Evolution for RELAMPAGO 2018 - 2019'])


%%
%figure
scatter(all_azimuth_correction_time, all_azimuth_correction_std)
datetick('x', 'HH:MM')
grid on 
%xlabel('HH:MM (UTC)')
datetick('x', 'mm/dd')
ylabel('Std. Azimuth Difference (deg.)')
title('Drift Evolution for RELAMPAGO 2018 - 2019')
