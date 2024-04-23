% DROPS for RELAMPAGO
% Ivan Arias
% 2019/05/10

% Includes fields to ingest data in pyart

input_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/DualRadar/Data/NetCDF_studyCase/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/Tallest_Storms/noDROPS/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_original/';
output_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/DualRadar/Data/DROPS/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/Tallest_Storms/DROPS/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/DROPS/';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/';
HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
sounding_file = ' ';
%config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_CHIVO_forCSAPR.v4.ini';
config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_CHIVO_forCHIVO.v4.ini';
%'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/Sounding_Bogota_20161113_00.txt';
%config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_CHIVO_forCSAPR.v4.ini';
%'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband.v1.ini';
%'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband_volumeMathing_v3.ini';
%'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband_volumeMathing_v2.ini';

membership_functions = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/membership_functions/';
directory = dir([input_path '/*col-radar*.nc']);
%directory = dir([input_path '/*csapr*.nc']);
%dir([input_path '/*col-radar*.nc']);
for i = 1:length(directory) 
    filename = directory(i).name;
    cmd = [HydroClass_exe ' ' input_path filename ' -o ' output_path filename ...
        ' -c ' config_file ' -s ' sounding_file ' -m VHS -d ' membership_functions];
    system(cmd);
    
    cmd = ['ncks -A -v fixed_angle ' input_path filename ' ' output_path filename];
    system(cmd);
    
    cmd = ['ncks -A -v sweep_mode ' input_path filename ' ' output_path filename];
    system(cmd);
    
    cmd = ['ncks -A -v time ' input_path filename ' ' output_path filename];
    system(cmd);
    
    
end