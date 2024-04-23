% 2021/11/10
% Zdr Calibration for CHIVO PRECIP
% Ivan Arias

addpath('/net/k2/storage/people/idariash/home/Utiles/PPI_Plotting');
addpath('/net/k2/storage/people/idariash/home/Utiles/Matlab');
addpath('/net/nasstore/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/code_matlab/utiles');
DataPath = '/net/k2/storage/people/idariash/home/Field_Campaigns/PRECIP/Zdr_calibration/NetCDF/20210731/';
Directory = dir([DataPath, '*.nc']);
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
    Z(Z < 0) = nan;
    Z(Z > 30) = nan;
%     ZDR(ZDR < -3) = nan;
%     ZDR(ZDR  > 3) = nan;
    
    
    Elv_range_indexing;
    
    Indexing_filter = Range < 50 & Elevation < 3 & Elevation > 0.9 & ...
        Azimuth > 90 & Azimuth < 180;
    
    Z_low_elev = Z(Indexing_filter)';
    
    ZDR_low_elev = ZDR(Indexing_filter)';
    
    Z_for_scatter = [Z_for_scatter Z_low_elev];
    ZDR_for_scatter = [ZDR_for_scatter ZDR_low_elev];
    disp(i)
 
end

Z= Z_for_scatter;
ZDR =  ZDR_for_scatter;