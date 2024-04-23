%filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/Case_210858Z/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';
%clear all
%close all
%filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/DualRadar/Data/DROPS/cfrad.20181214_020041.798_to_20181214_020556.800_col-radar_REL_PFAR360_SUR.nc';
filename = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2018/12/14/chivo.1b.20181214_020041.REL_PFAR360.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20181214/data/DROPS/corcsapr2cfrppiM1.a1.20181214.020004.nc';
azimuth_offset = 0;%7.4;

DBZ = ncread(filename, 'corrected_reflectivity');
ZDR = ncread(filename, 'corrected_differential_reflectivity');
RHOHV =  ncread(filename, 'corrected_copol_correlation_ratio');
KDP = ncread(filename, 'corrected_specific_differential_phase');
range = double(ncread(filename, 'range'))'/1e3;
azimuth = ncread(filename, 'azimuth') - azimuth_offset;
elevation = ncread(filename, 'elevation');



%% 
%Xo = 16;  % distance from the radar to get the vertical profile in km
delta = 0.47;
AZo = 156.5;%235;%335;%19;%257; %45;
Elv = 3.3;%2.5;

Index = abs(elevation-Elv) < 0.2 & abs(azimuth - AZo) < delta;
dbz = DBZ(:,Index);
zdr = ZDR(:,Index);
rhohv = RHOHV(:,Index);
kdp = KDP(:,Index);
x = range.*cosd(Elv);

zdr(dbz < 20) = nan;
rhohv(dbz < 20) = nan;
kdp(dbz < 20) = nan;
dbz(dbz < 20) = nan;


% X = Range.*cos(Elevation*pi/180);
% Y = Range.*sin(Elevation*pi/180);
% 
% Index = abs(Elevation-Elv) < 0.1 & abs(Azimuth - AZo) < 0.5;
% 
% y = Y(Index);
% x = X(Index);
% range = Range(Index);
% dbz = DBZ(Index);
% zdr = ZDR(Index);
% rhohv = RHOHV(Index);
% kdp = KDP(Index);
% 
% [x, I] = sort(x);
% dbz = dbz(I);
% zdr = zdr(I);
% rhohv = rhohv(I);
% kdp = kdp(I);
% 
% % zdr(zdr < -100) = nan;
% % plot(range, zdr)
% % %zdr(dbz < 10) = nan;
% % plot(range, zdr)
% % xlim([0, 56])
% 
% L = length(x) - mod(length(x),3);
% 
% x = x(1:L);
% dbz = dbz(1:L);
% zdr = zdr(1:L);
% rhohv = rhohv(1:L);
% kdp = kdp (1:L);
% 
% x = reshape(x,3, L/3);
% dbz = reshape(dbz,3, L/3);
% zdr = reshape(zdr,3, L/3);
% rhohv = reshape(rhohv,3, L/3);
% kdp = reshape(kdp,3, L/3);
% 
% x = nanmean(x);
% dbz = nanmean(dbz);
% zdr = nanmean(zdr);
% rhohv = nanmean(rhohv);
% kdp = nanmean(kdp);
% 
% dbz(dbz < -1000) = nan;
% zdr(zdr < -1000) = nan;
% rhohv(rhohv < -1000) = nan;
% kdp(kdp < -1000) = nan;
% 
% %filters
% zdr(dbz < 10) = nan;
% rhohv(dbz < 10) = nan;
% kdp(dbz < 10) = nan;
% dbz(dbz < 10) = nan;

% dbz(rhohv < 0.6) = nan;
% zdr(rhohv < 0.6) = nan;
% kdp(rhohv < 0.6) = nan;

%zdr = zdr + 0.7;

figure(2)
%title(['Vertical Profile at ' num2str(Xo) ' km'])

xx = -100:100;
yy_1 = 0*ones(1,length(xx));
yy_o = 0*ones(1,length(xx));
yy_2 = 0*ones(1,length(xx));

subplot(2,2,1)
plot(x, dbz)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
ylim([0,70]);
xlim([0,80]);
ylabel('(dBZ)')
xlabel('Distant from radar (km)')
grid on
title('Reflectivity')

subplot(2,2,2)
plot(x, zdr)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
ylim([-8, 8]);
xlim([0,80]);
ylabel('(dB)')
xlabel('Distant from radar (km)')
grid on
title('Differential Reflectivity')

subplot(2,2,3)
plot(x, rhohv)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
ylim([0.7, 1]);
xlim([0,80]);
ylabel('\rho_{hv} ')
xlabel('Distant from radar (km)')
grid on
title('Co-Polar Correlation')

subplot(2,2,4)
plot(x, kdp)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
ylim([-1, 8]);
xlim([0,80]);
ylabel('(deg/km)')
xlabel('Distant from radar (km)')
grid on
title('Specific Differential Phase')




