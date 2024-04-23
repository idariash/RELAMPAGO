% Ivan Arias
% 2020/12/17
clear all
addpath /net/denali/storage2/radar2/people/idariash/home/Utiles/PPI_Plotting
DataPath = '/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/rsc/out/preveous_1b';
'/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/rsc/out/1b.2_elbert';
'/net/denali/storage2/radar2/people/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/out/whitney_1b';



filename = '/net/denali/storage2/radar2/people/idariash/home/CSU/RELAMPAGO/analysis/WRF_Radar/out/chivo_grid_20181214T0300.nc';

directory = dir([DataPath '/*nc']);
L = length(directory);
ref_max_overTime = nan(40,L);
echoTop = nan(L,1);
%Time = nan(L,1);
smoth_filter =  ones(20,20,2); ones(2,2,2);
smoth_filter = 1/sum(sum(sum(smoth_filter)))*smoth_filter;
% smoth_filter2D = 1/4^2*ones(4,4);
% smoth_filter(:,:,2) = smoth_filter2D;
for i = 1:L
    filename = [DataPath '/' directory(i).name];
    ref = ncread(filename, 'corrected_reflectivity');
    x = ncread(filename, 'x')/1000;
    y = ncread(filename, 'y')/1000;
    z = ncread(filename, 'z')/1000;
    time = datetime(ncreadatt(filename,'/','start_datetime'), ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:SS''Z');
    
    [X,Y, Z] = meshgrid(x,y,z);
    R = sqrt(X.^2 + Y.^2);
    ref(R > 75) = nan;
    ref = convn(ref, smoth_filter, 'same');
    %ref(ref < 10) = nan;
    ref_max_overZ = max(max(ref));
    ref_max_overZ = reshape(ref_max_overZ, 40, 1);
    ref_max_overTime(:,i) = ref_max_overZ;
    echoTop(i) = max(z(ref_max_overZ > 18));
    Time(i) = time;
end

z = ncread(filename, 'z')/1000;

%ref_max_overTime = ref_max_overTime - 18;
ref_max_overTime(ref_max_overTime<-10) = nan;
time = 1:40;

figure
pcolor(Time,z,ref_max_overTime)
shading flat
load /net/denali/storage2/radar2/people/idariash/home/radartoolbox/colormaps/nws.mat
colormap(nwsZ)
caxis([-35, 75]);
h = colorbar;
ylabel(h, 'dBZ')

% hold on
% plot(Time, echoTop)
% hold off