addpath /net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting
filename_GPM = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/2A.GPM.DPR.V8-20180723.20190131-S222057-E235331.027991.V06A.HDF5';
filename_GR = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/data/NetCDF/20190131/corcsapr2cfrppiM1.a1.20190131.223003.nc';

figure

swath_number = 37;
%Read data
%ZcKu = double(h5read(filename_GPM,'/NS/PRE/zFactorMeasured'));
ZcKu = double(h5read(filename_GPM,'/NS/SLV/zFactorCorrected'));
ZcKu(ZcKu<15) = nan;
lat_ku = double(h5read(filename_GPM,'/NS/Latitude'));
lon_ku = double(h5read(filename_GPM,'/NS/Longitude')); 

ZcKa = double(h5read(filename_GPM,'/MS/SLV/zFactorCorrected'));
ZcKa(ZcKa<15) = nan;
lat_ka = double(h5read(filename_GPM,'/MS/Latitude'));
lon_ka = double(h5read(filename_GPM,'/MS/Longitude')); 


%select height
h_idx = 160;%162; %160~2Km from ground
sweep = 3;

%select scan number
scan_st =1;
scan_ed =7900;

%Plot GPM
figure(); hold on;
GPM = pcolor(lon_ku(:,scan_st:scan_ed),lat_ku(:,scan_st:scan_ed),squeeze(ZcKu(h_idx,:,scan_st:scan_ed)));
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
Z = ncread(filename_GR, 'reflectivity');
range = ncread(filename_GR, 'range');
elv = ncread(filename_GR, 'elevation');
azi = ncread(filename_GR, 'azimuth');
lat1 = ncread(filename_GR, 'latitude'); 
lon1 = ncread(filename_GR, 'longitude'); 

%azi = azi - 4; % correction of azimuth


elv = elv(1081:1440); %elv((sweep - 1)*360 + 1:360*sweep);
azi = azi(1081:1440); %azi((sweep - 1)*360 + 1:360*sweep);
Z_1st_sweep = Z(:,1081:1440); %Z((sweep - 1)*360*length(range) + 1: sweep*360*length(range)); %Z(1:360*length(range));
%Z_1st_sweep = vec2mat(Z_1st_sweep, length(range));
%Z_1st_sweep = Z_1st_sweep';
Z_1st_sweep(Z_1st_sweep < 10) = nan;
%Z_1st_sweep = Z_1st_sweep(range < 100e3, :);
%range = range(range<100e3);

e = referenceEllipsoid('WGS84');
[AZ,R]=meshgrid(azi,range);

[latout,lonout] = reckon(lat1,lon1,R,AZ,e);


refl = Z_1st_sweep;
HyC = pcolor(lonout,latout,double(refl));
units = 'HyC';
shading flat
%axis equal
color_map = 'jet'; 
title('2019/01/13 04:01 UTC','FontSize',10);
hc=colorbar;
xlabel(hc,'Reflectivity (dBZ)', 'FontSize',18);

HyC.FaceAlpha = 0.7;




lat =latout;
lon = lonout;

plot(lon1,lat1,'k*','MarkerSize',10,'LineWidth',5); 

plot(lonout(990,:),latout(990,:),'k','LineWidth',1);

plot(lonout(660,:),latout(660,:),'k','LineWidth',1);

%plot(lonout(500,:),latout(500,:),'k','LineWidth',1);

%plot(lonout(495,:),latout(495,:),'k','LineWidth',1);

%plot(lonout(300,:),latout(300,:),'k','LineWidth',1);

plot(lonout(330,:),latout(330,:),'k','LineWidth',1);

%plot(lonout(165,:),latout(165,:),'k','LineWidth',1);


xlabel('Longitude','FontSize',18,'FontWeight','bold');
set(gca,'fontsize',16,'FontWeight','bold');
ylabel('Latitude','FontSize',18,'FontWeight','bold');
set(gca,'fontsize',16,'FontWeight','bold');


xlim([-67, -62])
ylim([-33.75, -30.75])

hold off

%% ---------------------------------------------------------------

lat_swath = lat_ku(swath_number,:);
lon_swath = lon_ku(swath_number,:);
Lat_lon_index = -33.5 < lat_swath & lat_swath < -31.5 & -70 < lon_swath & lon_swath < -60;
ZcKu_RHI = ZcKu(65:176, swath_number, Lat_lon_index);
[r,~,c] = size(ZcKu_RHI);
ZcKu_RHI = reshape(ZcKu_RHI, r, c);
index_vcut = find(Lat_lon_index);
lat_rhi = lat_swath(Lat_lon_index);
lon_rhi = lon_swath(Lat_lon_index);
z = (r:-1:1)*0.125;
figure

pcolor(index_vcut,z,squeeze(ZcKu_RHI));
hold on 
%plot(lat1*ones(1,r),z,'k--','linewidth',1.5);
hold off

%figure

%pcolor(lon_rhi,z,squeeze(ZcKu_RHI));
colormap jet
shading flat
grid on
xlabel('Along track (scan index)')
ylabel('Height (km)')
title('DPR-Ku Vertical cut along bin 37')
colorbar
caxis([10 55]);
hc=colorbar;
xlabel(hc,'Corrected Reflectivity Ku (dBZ)');
%set(gca, 'XDir','reverse')
%xlim([-33, -31.5])

%% -------------------------------------------------


lat_swath = lat_ka(swath_number - 12,:);
lon_swath = lon_ka(swath_number - 12,:);
%Lat_lon_index = -33 < lat_swath & lat_swath < -31.5 & -65 < lon_swath & lon_swath < -63;
ZcKa_RHI = ZcKa(65:176, swath_number - 12, Lat_lon_index);
[r,~,c] = size(ZcKa_RHI);
ZcKa_RHI = reshape(ZcKa_RHI, r, c);
lat_rhi = lat_swath(Lat_lon_index);
lon_rhi = lon_swath(Lat_lon_index);
z = (r:-1:1)*0.125;
figure

pcolor(index_vcut,z,squeeze(ZcKa_RHI));
hold on 
%plot(lat1*ones(1,r),z,'k--','linewidth',1.5);
hold off

%figure

%pcolor(lon_rhi,z,squeeze(ZcKu_RHI));
colormap jet
shading flat
grid on
xlabel('Along track (scan index)')
ylabel('Height (km)')
title('DPR-Ka Vertical cut along bin 37')
colorbar
caxis([10 55]);
hc=colorbar;
xlabel(hc,'Corrected Reflectivity Ka (dBZ)');
%set(gca, 'XDir','reverse')
%xlim([-33, -31.5])

%------------------------------------------------------------

diff_freq_ratio = ZcKu_RHI - ZcKa_RHI;

figure

pcolor(index_vcut,z,squeeze(diff_freq_ratio));
hold on 
%plot(lat1*ones(1,r),z,'k--','linewidth',1.5);
hold off

%figure

%pcolor(lon_rhi,z,squeeze(ZcKu_RHI));
colormap jet
shading flat
grid on
xlabel('Along track (scan index)')
ylabel('Height (km)')
title('DPR-dfr Vertical cut along bin 37')
colorbar
caxis([-1 10]);
hc=colorbar;
xlabel(hc,'Diff. Freq. Ratio (dB)');
%set(gca, 'XDir','reverse')
%xlim([-33, -31.5])

