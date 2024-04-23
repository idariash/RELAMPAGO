% DROPS for RELAMPAGO
% Ivan Arias
% 2019/05/10

% Includes fields to ingest data in pyart

input_path = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/NetCDF_CHIVO/';
output_path = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/DROPS/CHIVO/';
HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
sounding_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/Sounding_Bogota_20161113_00.txt';
config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CHIVO_Cband.v3.ini';
membership_functions = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/membership_functions/';
directory = dir([input_path '/*.nc']);
for i = 1:length(directory) 
    filename = directory(i).name;
    cmd = [HydroClass_exe ' ' input_path filename ' -o ' output_path filename ...
        ' -c ' config_file ' -s ' sounding_file ' -m VHS -d ' membership_functions];
    system(cmd);
    
    cmd = ['ncks -A -v fixed_angle ' input_path filename ' ' output_path filename];
    system(cmd);
    
    
end
