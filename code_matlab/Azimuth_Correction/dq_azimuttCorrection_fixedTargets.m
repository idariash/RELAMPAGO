% Data Quality, Azimuth Error, fix target test
% 2019/07/24 Ivan based on Francesc suggestion
% 2019/07/29 find the fixed points

%% Read data
%data_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_HYDRO/20181114'; 
data_path = '/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181129';
directory = dir([data_path '/*HYDRO*']);
mask = ones(124, 360); %ones(166,360);
%mask = ones(166,360);
for i =1:6
    filename = [data_path '/' directory(i).name];
    Range = ncread(filename, 'range')/1e3; % range in km
    azimuth = ncread(filename, 'azimuth');
    elevation = ncread(filename, 'elevation');
    DBZ_TOT = ncread(filename, 'DBZ_TOT');
    VEL = ncread(filename, 'VEL');
    RHOHV = ncread(filename, 'RHOHV');
    

    % Filter
    figure 
    range_limit = 25;
    dbz_threshold = 45;
    rhohv_threshold = 0.6;
    %DBZ_TOT(RHOHV > rhohv_threshold) = nan; 
    dbz_tot = DBZ_TOT(Range < range_limit, elevation < 1);
    dbz_tot(dbz_tot < dbz_threshold) = nan;
    range = double(Range(Range < range_limit));
    azimuth = double(azimuth(elevation < 1));
    [azimuth, I] = sort(azimuth);
    dbz_tot = dbz_tot(:,I);
    mask_i = dbz_tot > dbz_threshold;
    mask = mask.*mask_i;
    %pcolor(azimuth, range, dbz_tot)
    imagesc(mask_i)
    shading flat
    grid on
%     colorbar
%     %colormap('jet')
%     xlabel('Azimuth (deg.)')
%     ylabel('Range (km)')
%     title('Total Reflectivity')
    
end
%%
imagesc(mask);
return
original_mask = mask;
figure
mask(mask < 0.5) = nan;
pcolor(azimuth, range, mask)
shading flat
grid on

%%

diff_mask = diff(original_mask, 1, 2);
figure
diff_mask(diff_mask > -0.5) = nan;
pcolor(diff_mask)
shading flat
grid on

return
filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_HYDRO/20181114/cfrad.20181114_140451.667_to_20181114_140541.892_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Hydro360/cfrad.20181110_235004.684_to_20181110_235055.038_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Hydro360/cfrad.20181110_234003.663_to_20181110_234054.825_col-radar_REL_HYDRO360_SUR.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Hydro360/cfrad.20181110_233004.529_to_20181110_233054.822_col-radar_REL_HYDRO360_SUR.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_Hydro360/cfrad.20181110_232005.136_to_20181110_232055.357_col-radar_REL_HYDRO360_SUR.nc';

Range = ncread(filename, 'range')/1e3; % range in km
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
DBZ_TOT = ncread(filename, 'DBZ_TOT');

%% Filter
figure 
range_limit = 25;
dbz_tot = DBZ_TOT(Range < range_limit, elevation < 1);
dbz_tot(dbz_tot < 50) = nan;
range = double(Range(Range < range_limit));
azimuth = double(azimuth(elevation < 1));
[azimuth, I] = sort(azimuth);
dbz_tot = dbz_tot(:,I);
pcolor(azimuth, range, dbz_tot)
shading flat
grid on
colorbar
colormap('jet')
xlabel('Azimuth (deg.)')
ylabel('Range (km)')
title('Total Reflectivity')
return

%% Isolate fixed points
dbz_tot_mask = dbz_tot > 40;
diff_dbz_mask = diff(dbz_tot_mask, 1, 2);
%diff_dbz_mask(diff_dbz_mask < 0.5) = nan;
figure
pcolor(diff_dbz_mask)
shading flat


