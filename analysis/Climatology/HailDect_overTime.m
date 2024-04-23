% Ivan Arias
% 2020/12/17
clear all

DataPath = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/mat_files';
addpath /net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/src
directory = dir([DataPath '/*.mat']);
L = length(directory);
HydroClass_max_overTime = nan(L,401);
%Time = nan(L,1);
smoth_filter = ones(1,1);
smoth_filter = 1/(sum(sum(smoth_filter)))*smoth_filter;
% smoth_filter2D = 1/4^2*ones(4,4);
% smoth_filter(:,:,2) = smoth_filter2D;
for i = 1:L
    filename = [DataPath '/' directory(i).name];
    load(filename);
%     HydroClass = ncread(filename, 'HydroClass');
%     x = ncread(filename, 'x')/1000;
%     y = ncread(filename, 'y')/1000;
%     z = ncread(filename, 'z')/1000;
    time = datetime(time, 'InputFormat', 'yyyyMMdd''T''HHmm');
    
    R = sqrt(x.^2 + y.^2);
    HydroClass(R > 100) = nan;
    HydroClass(HydroClass < 10) = nan;
    HydroClass(HydroClass > 11) = nan;
    HydroClass(max_ref<40) = nan;
    %HydroClass(HydroClass < 1) = 0;
    HydroClass = convn(HydroClass, smoth_filter, 'same');
    HydroClass_overX = max(HydroClass, [], 1);
    %HydroClass_overX = reshape(HydroClass_overX, 401, 1);
    HydroClass_max_overTime(i,:) = HydroClass_overX;
    Time(i) = time;
end

x = max(y, [], 2);

%HydroClass_max_overTime = HydroClass_max_overTime - 18;
%HydroClass_max_overTime(HydroClass_max_overTime<1) = nan;
figure
pcolor(x,Time,HydroClass_max_overTime)
% set(gca,'ColorScale','log')
shading flat
xlabel('South-North Distance to the CHIVO (km)')
ylabel('Time (UTC)')
xlim([-100 100])
hc=colorbar;
xlabel(hc,'Hail detection');
caxis([1 10]);
colormap(cmocean('rain'))
% russ_rain_colormap;
% %colormap(parula)
% %colormap(cmocean('rain'))
% colormap(cmap_data)


% hold on
% terrain
% hold off
% set(gca,'ColorScale','log')
