% Ivan Arias
% 2020/12/17
clear all
DataPath = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Climatology/mat_files/*.mat';
filename = '/net/denali/storage2/radar2/people/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/out/chivo_grid_20181214T0300.nc';

directory = dir(DataPath);
L = length(directory);
RainRate_max_overTime = nan(L,401);
echoTop = nan(L,1);
%Time = nan(L,1);
smoth_filter = ones(10,10,2);
smoth_filter = 1/sum(sum(sum(smoth_filter)))*smoth_filter;
% smoth_filter2D = 1/4^2*ones(4,4);
% smoth_filter(:,:,2) = smoth_filter2D;
for i = 1:L
    filename = ['./out/' directory(i).name];
    load(filename);
    RainRate = ncread(filename, 'RainRate');
    x = ncread(filename, 'x')/1000;
    y = ncread(filename, 'y')/1000;
    z = ncread(filename, 'z')/1000;
    time = datetime(ncreadatt(filename,'/','start_datetime'), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:SS''Z');
    
    [X,Y, Z] = meshgrid(x,y,z);
    R = sqrt(X.^2 + Y.^2);
    RainRate(R > 100) = nan;
    RainRate(RainRate < 0) = 0;
    RainRate = convn(RainRate, smoth_filter, 'same');
    RainRate(RainRate < 0) = nan;
    RainRate_overX = max(max(RainRate, [], 3));
    %RainRate_overX = reshape(RainRate_overX, 401, 1);
    RainRate_max_overTime(i,:) = RainRate_overX;
    Time(i) = time;
end

x = ncread(filename, 'x')/1000;

%RainRate_max_overTime = RainRate_max_overTime - 18;
RainRate_max_overTime(RainRate_max_overTime<1) = nan;
figure
pcolor(x,Time,RainRate_max_overTime)
set(gca,'ColorScale','log')
shading flat
russ_rain_colormap;
%colormap(parula)
%colormap(cmocean('rain'))
colormap(cmap_data)
caxis([1, 1000]);
h = colorbar;
ylabel(h, 'mm/h')

hold on
terrain
hold off
set(gca,'ColorScale','log')
