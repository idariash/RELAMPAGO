% Ivan Arias
% Plot Nexrad for San Juan Radar
% May-11/2018
% Colorado State University/Radar and Communcation Group

addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting');
addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab');
DataPath = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Calibration/Data/NetCDF/20190130/';
Directory = dir([DataPath, '*_P*']);
rhoHV_threshold = 0.9;
Z_for_scatter = [];
ZDR_for_scatter = [];
RadarType = 'Vaisala';
Campaign = 'RELAMPAGO';
disp(length(Directory))
for i = 1:length(Directory)
    close all
    filename = [DataPath, Directory(i).name];

    % get the fields
    Z = ncread(filename, 'DBZ');
    ZDR = ncread(filename, 'ZDR');
    %KDP = ncread(filename, 'KDP');
    RHOHV =  ncread(filename, 'RHOHV');
    range = double(ncread(filename, 'range'))'/1e3;
    azimuth = ncread(filename, 'azimuth');
    elevation = ncread(filename, 'elevation');
    ray_n_gates = ncread(filename, 'ray_n_gates');
    
    
    %RHOHV filter
    
    Z(RHOHV < 0.95) = nan;
    ZDR(RHOHV < 0.95) = nan;
    
    Elv_range_indexing;
    
    Z_low_elev = Z(Range < 50 & Elevation < 2 & Elevation > 1)';
    
    ZDR_low_elev = ZDR(Range < 50 & Elevation < 2 & Elevation > 1)';
    
    Z_for_scatter = [Z_for_scatter Z_low_elev];
    ZDR_for_scatter = [ZDR_for_scatter ZDR_low_elev];
    disp(i)
 
end

Z= Z_for_scatter;
ZDR =  ZDR_for_scatter;