% DROPS for RELAMPAGO
% Ivan Arias
% 2019/05/10

% Includes fields to ingest data in pyart

day = [18 19 20 21 27 28 02 03 06 07 08 09 10 13 14 15 17 18 22 23 24 25 26 28 29 30 31];
month = [12 12 12 12 12 12 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01 01];
year = [2018 2018 2018 2018 2018 2018 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019 2019];

for j = 1:27   
    input_path = ['/net/denali/storage/radar/RELAMPAGO/NetCDF/' num2str(year(j)) num2str(month(j), '%02.f') num2str(day(j), '%02.f') '/'];
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/';
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/'; 
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_original/';
    output_path = ['/net/denali/storage/radar/RELAMPAGO/DROPS/' num2str(year(j)) '/' num2str(month(j), '%02.f') '/' num2str(day(j), '%02.f') '/'];
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/';
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_1930/DROPS/';
    %'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/20190113/data/NetCDF_drops/';
    HydroClass_exe = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/hydroclass';
    sounding_file = ' ';
    %'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/Sounding_Bogota_20161113_00.txt';
    config_file = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_CHIVO_forCHIVO.v4.ini';
    %'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband.v1.ini';
    %'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband_volumeMathing_v3.ini';
    %'/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/RE_CSAPR_Cband_volumeMathing_v2.ini';

    membership_functions = '/net/denali/storage2/radar2/tmp/Ivan/Utiles/DROPS2/membership_functions/';
    directory = dir([input_path '/*RHI*.nc']);
    %dir([input_path '/*csapr*.nc']);
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
        return


    end
end