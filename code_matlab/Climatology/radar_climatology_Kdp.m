% Radar Climatology Statitistics
% Ivan Arias
% 2019/10/30

clear all
filename = '/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/analysis/relampago_statistics_Drops_Kdp.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
max_Kdp = T.max_Kdp;

%index_toFilter = echo_top > 2 & max_reflectivity > 0;
% echo_top(index_toFilter) = nan;
% max_reflectivity(index_toFilter) = nan;
% max_Kdp(index_toFilter) = nan;
% max_Zdr(index_toFilter) = nan;
% min_Zdr(index_toFilter) = nan;

start_date = datetime(2018,11,10);
end_date = datetime(2019,01, 31);
k = 1;
%min_Zdr(min_Zdr < -10) = nan;
for i = start_date:end_date    
    
    if length(max_Kdp(i < time_UTC & time_UTC < i + 1)) == 0
        continue
    end
    days(k) = i;
    maxKdp_day = max_Kdp(i < time_UTC & time_UTC < i + 1);
    maxKdp_mean(k) = nanmean(maxKdp_day);
    maxKdp_std(k) = nanstd(maxKdp_day);
    maxKdp_max(k) =  max(maxKdp_day);
    
    
    %echoTop_days(k,1:L) = echoTop_day; 
    k = k + 1;    
end
days = datenum(days);
return

figure
errorbar(days, echoTop_mean, echoTop_std);
datetick('x', 'dd')
hold on 
scatter(days, echoTop_max)
ylim([0 20])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Echo top height (km)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

figure
errorbar(days, maxRef_mean, maxRef_std);
datetick('x', 'dd')
hold on 
scatter(days, maxRef_max)
ylim([0 80])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Max. Reflectivity (dBZ)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

figure
errorbar(days, maxKdp_mean, maxKdp_std);
datetick('x', 'dd')
hold on 
scatter(days, maxKdp_max)
ylim([-3 17])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Max. K_{dp} (km)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

figure
errorbar(days, maxZdr_mean, maxZdr_std);
datetick('x', 'dd')
hold on 
scatter(days, maxZdr_max)
ylim([-3 10])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Max. Z_{dr} (dB)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

figure
errorbar(days, minZdr_mean, minZdr_std);
datetick('x', 'dd')
hold on 
scatter(days, minZdr_min)
%ylim([1 -10])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Min. Z_{dr} (dB)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

%boxplot(echoTop_days', days)
%datetick('x', 'dd-mmm')
