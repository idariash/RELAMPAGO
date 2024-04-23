% Ivan Arias
% 2019/12/20
% Get the value from CSAPR of a CHIVO measurement


r_chivo = 42;%37;
azi_chivo = 257;%234;
elv_chivo = 3.3;%2.3;

[r_csapr, azi_csapr, elv_csapr] = compute_csapr_coordinates(r_chivo, azi_chivo, elv_chivo);
%%
filename = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20181214/data/DROPS/corcsapr2cfrppiM1.a1.20181214.020004.nc';

DBZ = ncread(filename, 'corrected_reflectivity');
ZDR = ncread(filename, 'corrected_differential_reflectivity');
RHOHV =  ncread(filename, 'corrected_cross_correlation_ratio');
KDP = ncread(filename, 'corrected_specific_differential_phase');
range = double(ncread(filename, 'range'))'/1e3;
azimuth = ncread(filename, 'azimuth') - azimuth_offset;
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
r = Range.*cos(Elevation*pi/180); % projection over the ground
AZo = azi_csapr;
Elv = 2.3;

Index = abs(Elevation-Elv) < 0.5 & abs(Azimuth - AZo) < 1 & abs(r - r_csapr) < 1;
dbz = DBZ(Index);
zdr = ZDR(Index);
rhohv = RHOHV(Index);
kdp = KDP(Index);

