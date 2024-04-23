% Radar Climatology Statitistics
% Ivan Arias
% 2019/10/30

clear all
filename = '/top/students/GRAD/ECE/idariash/home/CSU//RELAMPAGO/analysis/Climatology/chivo_rainRate_hail.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
rain_accumulation = T.rain_accumulation;


index_toFilter = rain_accumulation > 3; %& max_reflectivity > 18 & max_reflectivity < 70;

start_date = datetime(2018,11,1);
end_date = datetime(2019,01, 31);
k = 1;
% min_Zdr(min_Zdr < -10) = nan;
for i = start_date:end_date
    rain_day = rain_accumulation(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    %maxRef_day = max_reflectivity(i < time_UTC & time_UTC < i + 1 & index_toFilter) + 1.7;
%     maxKdp_day = max_Kdp(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     maxZdr_day = max_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     minZdr_day = min_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    L = length(rain_day);    
    if L < 10
        continue
    end
    days(k) = i;
    rain_day_accumulation(k) = nansum(rain_day);
    
    k = k + 1;    
end


days = datenum(days);
rain_day_accumulation = rain_day_accumulation/(30000 * 6); % number of point in the grid time the hours
%figure
%scatter(days, echoTop_mean);

hold on 
scatter(days, rain_day_accumulation)
datetick('x', 'dd')
ylim([0 100])
xlim([datenum(start_date) datenum(end_date)])
ylabel('Rain Accumulation per day (mm)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
grid on
%-----------------------------------

clear all
filename = '/top/students/GRAD/ECE/idariash/home/CSU//RELAMPAGO/analysis/Climatology/csapr_rainRate_hail.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
rain_accumulation = T.rain_accumulation;


index_toFilter = rain_accumulation > 3; %& max_reflectivity > 18 & max_reflectivity < 70;

start_date = datetime(2018,11,1);
end_date = datetime(2019,01, 31);
k = 1;
% min_Zdr(min_Zdr < -10) = nan;
for i = start_date:end_date
    rain_day = rain_accumulation(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    %maxRef_day = max_reflectivity(i < time_UTC & time_UTC < i + 1 & index_toFilter) + 1.7;
%     maxKdp_day = max_Kdp(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     maxZdr_day = max_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     minZdr_day = min_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    L = length(rain_day);    
    if L < 10
        continue
    end
    days(k) = i;
    rain_day_accumulation(k) = nansum(rain_day);
    
    k = k + 1;    
end


days = datenum(days);
rain_day_accumulation = rain_day_accumulation/(30000 * 4); % number of point in the grid time the hours
%figure
%scatter(days, echoTop_mean);

hold on 
scatter(days, rain_day_accumulation)
datetick('x', 'dd')
ylim([0 100])
xlim([datenum(start_date) datenum(end_date)])
ylabel('Rain Accumulation per day (mm)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
grid on


%-------------------------------------
return
clear all
filename = '/top/students/GRAD/ECE/idariash/home/CSU//RELAMPAGO/analysis/Climatology/csapr_echoTop_distribution_unbias.xlsx';

'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_Drops_ppi_rhi.xlsx';

%'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_Drops_ppi_rhi.xlsx';
%'/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_all_Drops.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
echo_top = T.echo_top;
%max_reflectivity = T.max_ref;
% max_Kdp = T.max_Kdp;
% max_Zdr = T.max_Zdr;
% min_Zdr = T.min_Zdr;

index_toFilter = echo_top > 3; %& max_reflectivity > 18 & max_reflectivity < 70;
% echo_top(index_toFilter) = nan;
% max_reflectivity(index_toFilter) = nan;
% max_Kdp(index_toFilter) = nan;
% max_Zdr(index_toFilter) = nan;
% min_Zdr(index_toFilter) = nan;

start_date = datetime(2018,11,01);
end_date = datetime(2019,01, 31);
k = 1;
% min_Zdr(min_Zdr < -10) = nan;
for i = start_date:end_date
    echoTop_day = echo_top(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    %maxRef_day = max_reflectivity(i < time_UTC & time_UTC < i + 1 & index_toFilter) + 1.7;
%     maxKdp_day = max_Kdp(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     maxZdr_day = max_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     minZdr_day = min_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    L = length(echoTop_day);    
    if L < 5
        continue
    end
    days(k) = i;
    echoTop_mean(k) = nanmean(echoTop_day);
    echoTop_std(k) = nanstd(echoTop_day);
    echoTop_max(k) =  max(echoTop_day);
    
%     maxRef_mean(k) = nanmean(maxRef_day);
%     maxRef_std(k) = nanstd(maxRef_day, 1);
%     maxRef_max(k) =  max(maxRef_day);
    
%     maxKdp_mean(k) = nanmean(maxKdp_day);
%     maxKdp_std(k) = nanstd(maxKdp_day);
%     maxKdp_max(k) =  max(maxKdp_day);
%     
%     maxZdr_mean(k) = nanmean(maxZdr_day);
%     maxZdr_std(k) = nanstd(maxZdr_day);
%     maxZdr_max(k) =  max(maxZdr_day);
%     
%     minZdr_mean(k) = nanmean(minZdr_day);
%     minZdr_std(k) = nanstd(minZdr_day);
%     minZdr_min(k) =  min(minZdr_day);
    
    %echoTop_days(k,1:L) = echoTop_day; 
    k = k + 1;    
end


days = datenum(days);

%scatter(days, echoTop_mean);
echoTop_max = echoTop_max + 1.141;

scatter(days, echoTop_max)
datetick('x', 'dd')
ylim([0 20])
xlim([datenum(start_date) datenum(end_date)])
ylabel('Echo top height (km)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
grid on


%----------------------------------

clear all

filename = '/top/students/GRAD/ECE/idariash/home/CSU//RELAMPAGO/analysis/Climatology/rma_echoTop_distribution_unbias.xlsx';

'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_Drops_ppi_rhi.xlsx';

%'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_Drops_ppi_rhi.xlsx';
%'/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_all_Drops.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
echo_top = T.echo_top;
%max_reflectivity = T.max_ref;
% max_Kdp = T.max_Kdp;
% max_Zdr = T.max_Zdr;
% min_Zdr = T.min_Zdr;

index_toFilter = echo_top > 3; %& max_reflectivity > 18 & max_reflectivity < 70;
% echo_top(index_toFilter) = nan;
% max_reflectivity(index_toFilter) = nan;
% max_Kdp(index_toFilter) = nan;
% max_Zdr(index_toFilter) = nan;
% min_Zdr(index_toFilter) = nan;

start_date = datetime(2018,11,01);
end_date = datetime(2019,01, 31);
k = 1;
% min_Zdr(min_Zdr < -10) = nan;
for i = start_date:end_date
    echoTop_day = echo_top(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    %maxRef_day = max_reflectivity(i < time_UTC & time_UTC < i + 1 & index_toFilter) + 1.7;
%     maxKdp_day = max_Kdp(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     maxZdr_day = max_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
%     minZdr_day = min_Zdr(i < time_UTC & time_UTC < i + 1 & index_toFilter);
    L = length(echoTop_day);    
    if L < 5
        continue
    end
    days(k) = i;
    echoTop_mean(k) = nanmean(echoTop_day);
    echoTop_std(k) = nanstd(echoTop_day);
    echoTop_max(k) =  max(echoTop_day);
    
%     maxRef_mean(k) = nanmean(maxRef_day);
%     maxRef_std(k) = nanstd(maxRef_day, 1);
%     maxRef_max(k) =  max(maxRef_day);
    
%     maxKdp_mean(k) = nanmean(maxKdp_day);
%     maxKdp_std(k) = nanstd(maxKdp_day);
%     maxKdp_max(k) =  max(maxKdp_day);
%     
%     maxZdr_mean(k) = nanmean(maxZdr_day);
%     maxZdr_std(k) = nanstd(maxZdr_day);
%     maxZdr_max(k) =  max(maxZdr_day);
%     
%     minZdr_mean(k) = nanmean(minZdr_day);
%     minZdr_std(k) = nanstd(minZdr_day);
%     minZdr_min(k) =  min(minZdr_day);
    
    %echoTop_days(k,1:L) = echoTop_day; 
    k = k + 1;    
end


days = datenum(days);

%scatter(days, echoTop_mean);

echoTop_max = echoTop_max + 0.484;

scatter(days, echoTop_max)
datetick('x', 'dd')
ylim([0 20])
xlim([datenum(start_date) datenum(end_date)])
ylabel('Echo top height (km)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
grid on


%-----------------------------------
return

figure
scatter(days, maxRef_mean);
datetick('x', 'dd')
hold on 
scatter(days, maxRef_max)
ylim([0 80])
xlim([datenum(start_date - 1) datenum(end_date + 1)])
ylabel('Max. Reflectivity (dBZ)')
legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
    'Location','northwest')
grid on

% figure
% errorbar(days, maxKdp_mean, maxKdp_std);
% datetick('x', 'dd')
% hold on 
% scatter(days, maxKdp_max)
% ylim([-3 17])
% xlim([datenum(start_date - 1) datenum(end_date + 1)])
% ylabel('Max. K_{dp} (km)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
% grid on
% 
% figure
% errorbar(days, maxZdr_mean, maxZdr_std);
% datetick('x', 'dd')
% hold on 
% scatter(days, maxZdr_max)
% ylim([-3 10])
% xlim([datenum(start_date - 1) datenum(end_date + 1)])
% ylabel('Max. Z_{dr} (dB)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
% grid on
% 
% figure
% errorbar(days, minZdr_mean, minZdr_std);
% datetick('x', 'dd')
% hold on 
% scatter(days, minZdr_min)
% %ylim([1 -10])
% xlim([datenum(start_date - 1) datenum(end_date + 1)])
% ylabel('Min. Z_{dr} (dB)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
% grid on

%boxplot(echoTop_days', days)
%datetick('x', 'dd-mmm')