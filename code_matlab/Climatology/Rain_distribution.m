% Radar Climatology Statitistics
% Ivan Arias
% 2019/10/30

clear all
filename = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_echoTop_distribution_100km_ppi.xlsx';
'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_echoTop_distribution_Dec14.xlsx';
'/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_Drops_ppi_rhi.xlsx';
%'/net/denali/storage2/radar2/tmp/Ivan/CSU/RELAMPAGO/analysis/Climatology/relampago_statistics_all_Drops.xlsx';

T = readtable(filename);

time_UTC = datetime(T.time_UTC, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
echo_top = T.echo_top;
x_echo_top = T.y_echo_top;

echo_top(echo_top < 4) = nan;

%max_reflectivity = T.max_ref;
% max_Kdp = T.max_Kdp;
% max_Zdr = T.max_Zdr;
% min_Zdr = T.min_Zdr;

index_toFilter = echo_top > 3;% & max_reflectivity > 18 & max_reflectivity < 70;
% echo_top(index_toFilter) = nan;
% max_reflectivity(index_toFilter) = nan;
% max_Kdp(index_toFilter) = nan;
% max_Zdr(index_toFilter) = nan;
% min_Zdr(index_toFilter) = nan;

start_date = datetime(2018,11,10);
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
    if L < 10
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

return

% Hovmoller 3 month

figure
start_date = datetime(2019,01, 20);
end_date = datetime(2019, 01, 31);
time_local = time_UTC;
scatter(x_echo_top,time_UTC,25, echo_top,'filled')
xlabel('North-South Distance to the CHIVO (km) ')
ylabel('Time (UTC) ')
xlim([-100 100])
ylim([start_date end_date])
colormap(parula(16))
hc=colorbar;
xlabel(hc,'Echo top height (km)');
caxis([4 20]);


% Diurnal cycle
figure
time_hour = mod(hour(time_UTC) - 3 + minute(time_UTC)./60, 24);
time_hour(time_hour < 12) = time_hour(time_hour < 12) + 24;
scatter(x_echo_top,time_hour,25, echo_top,'filled')
xlabel('North-South Distance to the CHIVO (km) ')
ylabel('Local time ')
xlim([-100 100])
ylim([12 36])
yticks([12 16 20 24 28 32 36])
yticklabels([12 16 20 0 4 8 12])
colormap(parula(10))
hc=colorbar;
xlabel(hc,'Echo top height (km)');
caxis([10 20]);


% Max echotop
figure
days = datenum(days);
% scatter(days, echoTop_mean);
hold on 
scatter(days, echoTop_max)
datetick('x', 'dd')
ylim([0 20])
xlim([datenum(start_date - 1) datenum(end_date +1)])
ylabel('Max. echo top height (km)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
%grid on



% 
% figure
% scatter(days, maxRef_mean);
% datetick('x', 'dd')
% hold on 
% scatter(days, maxRef_max)
% ylim([0 80])
% xlim([datenum(start_date - 1) datenum(end_date + 1)])
% ylabel('Max. Reflectivity (dBZ)')
% legend({'Distribution of the max. value observed every 10 min', 'Max. value observed during the whole day'},...
%     'Location','northwest')
% grid on

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