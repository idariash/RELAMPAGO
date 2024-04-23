
filename = '/net/denali/storage/radar/RELAMPAGO/NetCDF/20190125/cfrad.20190125_213004.545_to_20190125_213024.545_col-radar_BIRDBATH_SUR.nc';
%filename = '/net/denali/storage/radar/RELAMPAGO/NetCDF/20190125/cfrad.20190125_213026.652_to_20190125_213046.652_col-radar_BIRDBATH_SUR.nc';
range = ncread(filename, 'range')/1e3;
azimuth = ncread(filename, 'azimuth');
DBZ = ncread(filename, 'DBZ');
VEL = ncread(filename, 'VEL');
dbz = nanmean(DBZ, 2);
vel = nanmean(VEL, 2);
figure
plot(dbz, range)
ylim([0, 20])

figure
plot(vel, range)
ylim([0, 10])
xlim([-15, 15])
xlabel('Vertical velocity (m/s)')
ylabel('Height (km)')


figure

plot(azimuth, VEL(111,:)) %41
ylim([-15, 15])
xlim([0, 360])
xlabel('Azimuth (deg)')
ylabel('Velocity (m/s)')
