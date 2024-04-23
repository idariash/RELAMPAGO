% Ivan Arias
% 2020/12/17
clear all

DataPath = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/mat_files';
addpath /net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/src
directory = dir([DataPath '/*.mat']);
L = length(directory);
RainRate_max_overTime = nan(L,401);
%Time = nan(L,1);
smoth_filter = ones(1,1);
smoth_filter = 1/(sum(sum(smoth_filter)))*smoth_filter;
% smoth_filter2D = 1/4^2*ones(4,4);
% smoth_filter(:,:,2) = smoth_filter2D;
for i = 1:L
    filename = [DataPath '/' directory(i).name];
    load(filename);
%     RainRate = ncread(filename, 'RainRate');
%     x = ncread(filename, 'x')/1000;
%     y = ncread(filename, 'y')/1000;
%     z = ncread(filename, 'z')/1000;
    time = datetime(time, 'InputFormat', 'yyyyMMdd''T''HHmm');
    
    R = sqrt(x.^2 + y.^2);
    RainRate(R > 100) = nan;
    RainRate(RainRate < 4) = nan;
    RainRate(RainRate > 100) = nan;
    %RainRate(RainRate < 1) = 0;
    RainRate = convn(RainRate, smoth_filter, 'same');
    RainRate_overX = max(RainRate, [], 1);
    %RainRate_overX = reshape(RainRate_overX, 401, 1);
    RainRate_max_overTime(i,:) = RainRate_overX;
    Time(i) = time;
end

x = max(y, [], 2);

%RainRate_max_overTime = RainRate_max_overTime - 18;
%RainRate_max_overTime(RainRate_max_overTime<1) = nan;
figure
pcolor(x,Time,RainRate_max_overTime)
% set(gca,'ColorScale','log')
shading flat
xlabel('South-North Distance to the CHIVO (km)')
ylabel('Time (UTC)')
xlim([-100 100])
hc=colorbar;
xlabel(hc,'Rain Rate (mm/h)');
caxis([1 100]);
colormap(cmocean('rain'))
% russ_rain_colormap;
% %colormap(parula)
% %colormap(cmocean('rain'))
% colormap(cmap_data)


% hold on
% terrain
% hold off
% set(gca,'ColorScale','log')
