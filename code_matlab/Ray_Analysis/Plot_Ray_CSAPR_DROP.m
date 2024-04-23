% Ivan Arias
% Mayo 10/2019
% Ray over the strom 

filename = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/NetCDF_CHIVO/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/NetCDF/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/corcsapr2cfrhsrhiM1.a1.20190125.193715.nc';
DBZ = ncread(filename, 'DBZ');
ZDR = ncread(filename, 'ZDR');
%uncorrected_reflectivity_h = ncread(filename, 'uncorrected_reflectivity_h');
range = ncread(filename, 'range')'/1e3;
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
ray_n_gates = ncread(filename, 'ray_n_gates');
delta = 0.5;
Azi =250;
Elv = 2; 
zdr_bias = 0.7;

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


Index = abs(Azimuth - Azi) < delta & abs(Elevation - Elv) < delta;

dbz = DBZ(Index);
zdr = ZDR(Index);
range = Range(Index);
[range, I] = sort(range);
dbz = dbz(I);
zdr = zdr(I);

%------------
L = length(range) - mod(length(range),10);

range = range(1:L);
dbz = dbz(1:L);
zdr = zdr(1:L);

range = reshape(range,10, L/10);
dbz = reshape(dbz,10, L/10);
zdr = reshape(zdr,10, L/10);

zdr = zdr - zdr_bias;

range = nanmean(range);
dbz = nanmean(dbz);
zdr = nanmean(zdr);

dbz(dbz < -100) = nan;
zdr(zdr < -100) = nan;

zdr(dbz < 0 ) = nan;
dbz(dbz < 0) = nan;
x_range = range*cos(Elv*pi/180);

%-------------
figure(1)
plot(x_range, dbz)
hold on 

figure(2)
plot(x_range, zdr)
hold on 


%----------------------------------

filename = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/DROPS/CHIVO/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/DROPS/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/corcsapr2cfrhsrhiM1.a1.20190125.193715.nc';
DBZ = ncread(filename, 'corrected_reflectivity');
ZDR = ncread(filename, 'corrected_differential_reflectivity');
%uncorrected_reflectivity_h = ncread(filename, 'uncorrected_reflectivity_h');
range = ncread(filename, 'range')'/1e3;
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
ray_n_gates = ncread(filename, 'ray_n_gates');
delta = 0.5;

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

Index = abs(Azimuth - Azi) < delta & abs(Elevation - Elv) < delta;

dbz = DBZ(Index);
zdr = ZDR(Index);
range = Range(Index);
[range, I] = sort(range);
dbz = dbz(I);
zdr = zdr(I);

%------------
L = length(range) - mod(length(range),10);

range = range(1:L);
dbz = dbz(1:L);
zdr = zdr(1:L);

range = reshape(range,10, L/10);
dbz = reshape(dbz,10, L/10);
zdr = reshape(zdr,10, L/10);

range = nanmean(range);
dbz = nanmean(dbz);
zdr = nanmean(zdr);


dbz(dbz < 0) = nan;
zdr(zdr < 0) = nan;
x_range = range*cos(Elv*pi/180);

%-------------
figure(1)
plot(x_range, dbz)
xlim([0, 60])
ylim([10, 60])
xlabel('Range Projection over X (km)')
ylabel('Reflectivity (dBZ)')
grid on
hold off 

figure(2)
plot(x_range, zdr)
xlim([0, 60])
ylim([0, 6])
xlabel('Range Projection over X (km)')
ylabel('Diff. Reflectivity (dB)')
grid on
hold off 
