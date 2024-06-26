% Ivan Arias
% Mayo 10/2019
% Ray over the strom 

filename = '/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181114/cfrad.20181114_135139.736_to_20181114_135816.409_col-radar_REL_PPIFAR36_SUR.nc';
%'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181122/cfrad.20181122_175004.474_to_20181122_175041.549_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181116/cfrad.20181116_213614.669_to_20181116_213810.628_col-radar_REL_SURV_SUR.nc';
%'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181114/cfrad.20181114_140451.667_to_20181114_140541.892_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181111/cfrad.20181111_172005.425_to_20181111_172055.656_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181110/cfrad.20181110_203015.927_to_20181110_203106.222_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20190120_095042.141_to_20190120_095556.289_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20190116_094045.163_to_20190116_094600.466_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20190109_094042.998_to_20190109_094558.333_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20190104_094044.800_to_20190104_094558.784_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20181228_093042.269_to_20181228_093557.095_col-radar_REL_PFAR360_SUR.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20190125_095048.467_to_20190125_095604.801_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20181211_230043.127_to_20181211_230559.451_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20181201_092042.348_to_20181201_092556.859_col-radar_REL_PFAR360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Solar/cfrad.20181125_225008.928_to_20181125_225547.068_col-radar_PPICOWH360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/NetCDF_CHIVO/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/NetCDF/cfrad.20190125_193422.355_to_20190125_193746.772_col-radar_REL_RHI45_RHI.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/corcsapr2cfrhsrhiM1.a1.20190125.193715.nc';
DBZ = ncread(filename, 'DBZ_TOT');
%ZDR = ncread(filename, 'ZDR');
%uncorrected_reflectivity_h = ncread(filename, 'uncorrected_reflectivity_h');
range = ncread(filename, 'range')'/1e3;
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
%ray_n_gates = ncread(filename, 'ray_n_gates');
delta = 0.49;

% Range = [];
% Elevation = [];
% Azimuth = [];
% for I = 1:length(ray_n_gates)
%     Range = [Range range(1:ray_n_gates(I))];
%     
% end
% 
% for I = 1:length(elevation)
%     elevation_n_gates = elevation(I)*ones(1,ray_n_gates(I));
%     azimuth_n_gates = azimuth(I)*ones(1,ray_n_gates(I));
%     Elevation = [Elevation elevation_n_gates];
%     Azimuth = [Azimuth azimuth_n_gates];
% end

%%
Azi = 352;
Elv = 0.5; 
%Index = abs(azimuth - Azi) < delta & abs(elevation - Elv) < delta;

dbz = DBZ(:, abs(azimuth - Azi) < delta);
dbz = nanmean(dbz,2)';
%zdr = ZDR(Index);
%range = Range(Index);
%[range, I] = sort(range);
%dbz = dbz(I);
%zdr = zdr(I);

%------------
L = length(range) - mod(length(range),10);

range = range(1:L);
dbz = dbz(1:L);
%zdr = zdr(1:L);

range = reshape(range,10, L/10);
dbz = reshape(dbz,10, L/10);
%zdr = reshape(zdr,10, L/10);

range = nanmean(range);
dbz = nanmean(dbz);
%zdr = nanmean(zdr);

dbz(dbz < -100) = nan;
%zdr(zdr < -100) = nan;

%zdr(dbz < 0 ) = nan;
dbz(dbz < 0) = nan;
%x_range = range*cos(Elv*pi/180);

%-------------
figure
plot(range, dbz)
xlabel('range (km)')
ylabel('Total Reflectivity (dBZ)')
ylim([0,45])
xlim([0,150])
title(['Azimuth: ' num2str(Azi) '° | Elevation: ' num2str(Elv) '°'])
grid on
%hold on 

return

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

zdr = zdr + 0.7; %Bias correction

dbz(dbz < -100) = nan;
zdr(zdr < -100) = nan;
x_range = range*cos(Elv*pi/180);

%-------------
figure(1)
plot(x_range, dbz)
xlim([20, 60])
ylim([10, 60])
xlabel('Range Projection over X (km)')
ylabel('Reflectivity (dBZ)')
grid on
hold off 

figure(2)
plot(x_range, zdr)
xlim([20, 60])
ylim([-1, 7])
xlabel('Range Projection over X (km)')
ylabel('Diff. Reflectivity (dB)')
grid on
hold off 
