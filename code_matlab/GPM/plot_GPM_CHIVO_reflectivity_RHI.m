addpath /net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting
filename_GPM = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM_CSAPR/GPM/2A.GPM.DPR.V8-20180723.20190113-S034700-E051932.027699.V06A.HDF5'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/2A.GPM.DPR.V8-20180723.20181206-S040258-E053532.027108.V06A.HDF5';
filename_GR = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM_CSAPR/GPM/DROPS/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20181206/cfrad.20181206_052042.951_to_20181206_052559.414_col-radar_REL_PFAR360_SUR_DROPS.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/cfrad.20181206_051842.000_to_20181206_052007.914_1_SUR.nc';

figure

swath_number = 13;
%Read data
Zmku = double(h5read(filename_GPM,'/NS/PRE/zFactorMeasured'));
Zmku(Zmku<15) = nan;
lat_ku = double(h5read(filename_GPM,'/NS/Latitude'));
lon_ku = double(h5read(filename_GPM,'/NS/Longitude')); 

%select height
h_idx = 145;%162; %160~2Km from ground
sweep = 3;

%select scan number
scan_st =1;
scan_ed =7900;

%Plot GPM
figure(); hold on;
GPM = pcolor(lon_ku(:,scan_st:scan_ed),lat_ku(:,scan_st:scan_ed),squeeze(Zmku(h_idx,:,scan_st:scan_ed)));
colorbar;
caxis([10 50]);
xlabel('Longitude','fontsize',12);
ylabel('Latitude','fontsize',12); 
plot(lon_ku(1,scan_st:scan_ed),lat_ku(1,scan_st:scan_ed),'k','linewidth',1.5);
plot(lon_ku(49,scan_st:scan_ed),lat_ku(49,scan_st:scan_ed),'k','linewidth',1.5);
%plot(lon_ku(13,scan_st:scan_ed),lat_ku(13,scan_st:scan_ed),'k--','linewidth',1.5);
plot(lon_ku(swath_number,scan_st:scan_ed),lat_ku(swath_number,scan_st:scan_ed),'k--','linewidth',1.5);
plot_google_map('APIKey','AIzaSyB9LuHNhTyshkiPqCUNeio903DCrI-8v5U'); %get this function from interent. Its available to public. 
shading flat
colormap jet

%Plot Radar
Z = ncread(filename_GR, 'corrected_reflectivity');
range = ncread(filename_GR, 'range');
elv = ncread(filename_GR, 'elevation');
azi = ncread(filename_GR, 'azimuth');
lat1 = ncread(filename_GR, 'latitude'); 
lon1 = ncread(filename_GR, 'longitude'); 



elv = elv((sweep - 1)*360 + 1:360*sweep);
azi = azi((sweep - 1)*360 + 1:360*sweep);
Z_1st_sweep = Z((sweep - 1)*360*length(range) + 1: sweep*360*length(range)); %Z(1:360*length(range));
Z_1st_sweep = vec2mat(Z_1st_sweep, length(range));
Z_1st_sweep = Z_1st_sweep';
Z_1st_sweep(Z_1st_sweep < 10) = nan;
%Z_1st_sweep = Z_1st_sweep(range < 100e3, :);
%range = range(range<100e3);

azi = azi - 14;

e = referenceEllipsoid('WGS84');
[AZ,R]=meshgrid(azi,range);

[latout,lonout] = reckon(lat1,lon1,R,AZ,e);


refl = Z_1st_sweep;
HyC = pcolor(lonout,latout,double(refl));
units = 'HyC';
shading flat
%axis equal
color_map = 'jet'; 
title('2018/12/06 05:22 UTC','FontSize',10);
hc=colorbar;
xlabel(hc,'Reflectivity (dBZ)', 'FontSize',18);

HyC.FaceAlpha = 0.7;




lat =latout;
lon = lonout;

plot(lon1,lat1,'k*','MarkerSize',10,'LineWidth',5); 

plot(lonout(660,:),latout(660,:),'k','LineWidth',1);

%plot(lonout(500,:),latout(500,:),'k','LineWidth',1);

plot(lonout(495,:),latout(495,:),'k','LineWidth',1);

%plot(lonout(300,:),latout(300,:),'k','LineWidth',1);

plot(lonout(330,:),latout(330,:),'k','LineWidth',1);

plot(lonout(165,:),latout(165,:),'k','LineWidth',1);


xlabel('Longitude','FontSize',18,'FontWeight','bold');
set(gca,'fontsize',16,'FontWeight','bold');
ylabel('Latitude','FontSize',18,'FontWeight','bold');
set(gca,'fontsize',16,'FontWeight','bold');


xlim([-66, -62])
ylim([-33, -30])

hold off

lat_swath = lat_ku(swath_number,:);
lon_swath = lon_ku(swath_number,:);
Lat_lon_index = -33 < lat_swath & lat_swath < -30 & -65 < lon_swath & lon_swath < -64;
ZmKu_RHI = Zmku(80:161, swath_number, Lat_lon_index);
[r,~,c] = size(ZmKu_RHI);
ZmKu_RHI = reshape(ZmKu_RHI, r, c);
lat_rhi = lat_swath(Lat_lon_index);
lon_rhi = lon_swath(Lat_lon_index);
z = (r:-1:1)*0.25;
figure

pcolor(lat_rhi,z,squeeze(ZmKu_RHI));
hold on 
plot(-31.63*ones(1,r),z,'k--','linewidth',1.5);
hold off

%figure

%pcolor(lon_rhi,z,squeeze(ZmKu_RHI));
colormap jet
shading flat
grid on
xlabel('Latitude Projection')
ylabel('Height (km)')
colorbar
hc=colorbar;
xlabel(hc,'Reflectivity (dBZ)');
