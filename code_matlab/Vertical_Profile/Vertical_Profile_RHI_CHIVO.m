filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/Case_210858Z/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';

%filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/PyArt_test/DROPS_write_cfradial_0/cfrad.20190125_210858.592_to_20190125_211009.717_col-radar_REL_RHI30_RHI.nc';

Xo = 16;  % distance from the radar to get the vertical profile in km
delta = 0.5;
AZo = 282.5;

DBZ = ncread(filename, 'DBZ');
ZDR = ncread(filename, 'ZDR');
RHOHV =  ncread(filename, 'RHOHV');
KDP = ncread(filename, 'KDP');
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


figure(1)
%title(['Vertical Profile at ' num2str(Xo) ' km'])

xx = -100:100;
yy_1 = 8*ones(1,length(xx));
yy_o = 9*ones(1,length(xx));
yy_2 = 10*ones(1,length(xx));

subplot(2,2,1)
hold on
plot(dbz, y)

plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([0,70]);
ylim([8,10]);
xlabel('(dBZ)')
ylabel('Hieght(Km)')
grid on
title('Reflectivity')

subplot(2,2,2)
hold on

plot(zdr, y)
hold on
plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([-2, 5]);
ylim([8,10]);
xlabel('(dB)')
ylabel('Hieght(Km)')
grid on
title('Differential Reflectivity')

subplot(2,2,3)
hold on
plot(rhohv, y)

plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([0.5, 1]);
ylim([8,10]);
xlabel('\rho_{hv} ')
ylabel('Hieght(Km)')
grid on
title('Cross Polar Correlation')

subplot(2,2,4)
hold on
plot(kdp, y)

plot(xx, yy_1,'--');
plot(xx, yy_o,'r');
plot(xx, yy_2,'--');
hold off
xlim([-1, 8]);
ylim([8,10]);
xlabel('(deg/km)')
ylabel('Hieght(Km)')
grid on
title('Specific Differential Phase')




