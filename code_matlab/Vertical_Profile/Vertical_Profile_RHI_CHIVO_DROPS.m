%filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/Case_210858Z/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';
%clear all
%close all
filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Tallest_Storms/DROPS/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Tallest_Storms/DROPS/cfrad.20181214_020602.120_to_20181214_020923.124_col-radar_REL_RHI45_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Tallest_Storms/DROPS/cfrad.20181110_213530.450_to_20181110_213713.101_col-radar_REL_RHI45_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/PyArt_test/DROPS_write_cfradial_0/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';



DBZ = ncread(filename, 'corrected_reflectivity');
ZDR = ncread(filename, 'corrected_differential_reflectivity');
RHOHV =  ncread(filename, 'corrected_cross_correlation_ratio');
KDP = ncread(filename, 'corrected_specific_differential_phase');
range = double(ncread(filename, 'range'))'/1e3;
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
ray_n_gates = ncread(filename, 'ray_n_gates');

Range = [];
Elevation = [];
Azimuth = [];
for I = 1:length(ray_n_gates)
    Range = [Range range(1:ray_n_gates(I))];
    
end

for I = 1:length(elevation)
    elevation_n_gates = elevation(I)*ones(1,ray_n_gates(I));
    azimuth_n_gates = azimuth(I)*ones(1,ray_n_gates(I));
    Elevation = [Elevation elevation_n_gates];
    Azimuth = [Azimuth azimuth_n_gates];
end

%% 
Xo = 16;  % distance from the radar to get the vertical profile in km
delta = 0.5;
AZo = 285;

X = Range.*cos(Elevation*pi/180);
Y = Range.*sin(Elevation*pi/180);

Index = abs(X-Xo) < delta & abs(Azimuth - AZo) < delta;

y = Y(Index);
x = X(Index);
dbz = DBZ(Index);
zdr = ZDR(Index);
rhohv = RHOHV(Index);
kdp = KDP(Index);

[y, I] = sort(y);
dbz = dbz(I);
zdr = zdr(I);
rhohv = rhohv(I);
kdp = kdp(I);

L = length(y) - mod(length(y),10);

y = y(1:L);
dbz = dbz(1:L);
zdr = zdr(1:L);
rhohv = rhohv(1:L);
kdp = kdp (1:L);

y = reshape(y,10, L/10);
dbz = reshape(dbz,10, L/10);
zdr = reshape(zdr,10, L/10);
rhohv = reshape(rhohv,10, L/10);
kdp = reshape(kdp,10, L/10);

y = nanmean(y);
dbz = nanmean(dbz);
zdr = nanmean(zdr);
rhohv = nanmean(rhohv);
kdp = nanmean(kdp);

dbz(dbz < -1000) = nan;
zdr(zdr < -1000) = nan;
rhohv(rhohv < -1000) = nan;
kdp(kdp < -1000) = nan;

%zdr = zdr + 0.7;

figure(1)
%title(['Vertical Profile at ' num2str(Xo) ' km'])

xx = -100:100;
yy_1 = 0*ones(1,length(xx));
yy_o = 0*ones(1,length(xx));
yy_2 = 0*ones(1,length(xx));

subplot(2,2,1)
plot(dbz, y)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([0,70]);
ylim([0,20]);
xlabel('(dBZ)')
ylabel('Hieght(Km)')
grid on
title('Reflectivity')

subplot(2,2,2)
plot(zdr, y)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([-2, 8]);
ylim([0,20]);
xlabel('(dB)')
ylabel('Hieght(Km)')
grid on
title('Differential Reflectivity')

subplot(2,2,3)
plot(rhohv, y)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([0.5, 1]);
ylim([0,20]);
xlabel('\rho_{hv} ')
ylabel('Hieght(Km)')
grid on
title('Cross Polar Correlation')

subplot(2,2,4)
plot(kdp, y)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([-1, 8]);
ylim([0,20]);
xlabel('(deg/km)')
ylabel('Hieght(Km)')
grid on
title('Specific Differential Phase')




