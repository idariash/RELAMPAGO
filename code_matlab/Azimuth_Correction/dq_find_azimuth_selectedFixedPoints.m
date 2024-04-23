% Data Quality, Find Azimuth of the fixed points
% 2019/07/24 Ivan based on Francesc suggestion
% 

data_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Rotation/NetCDF_HYDRO/20181114';
directory = dir([data_path '/*HYDRO*']);

%% fixed Points were found by visutal inspection
fixedPoint_azimuthIndex = [285, 191, 239, 227, 296];
fixedPoint_rangeIndex = [122, 108, 99, 42, 11];
fixedPoint_azimuth = zeros(length(directory), 5);
fixedPoint_range = zeros(length(directory), 5);
fixedPoint_dbz = zeros(length(directory), 5);
fixedPoint_vel = zeros(length(directory), 5);
fixedPoint_rhohv = zeros(length(directory), 5);

for i =1:length(directory)
    filename = [data_path '/' directory(i).name];
    range = ncread(filename, 'range')/1e3; % range in km
    azimuth = ncread(filename, 'azimuth');
    elevation = ncread(filename, 'elevation');
    DBZ_TOT = ncread(filename, 'DBZ_TOT');
    VEL = ncread(filename, 'VEL');
    RHOHV = ncread(filename, 'RHOHV');

    % Sort azimuth and moments
    azimuth = azimuth(elevation < 1); % First sweep of HYDRO
    dbz_tot = DBZ_TOT(:, elevation < 1);
    vel = VEL(:, elevation < 1);
    rhohv = RHOHV(:, elevation < 1);
    
    [azimuth, I] = sort(azimuth);
    dbz_tot = dbz_tot(:,I);
    vel = vel(:, I);
    rhohv = rhohv(:,I);
    
    for j = 1:length(fixedPoint_azimuthIndex) % to average
        
        fixedPoint_azimuth(i,j) = azimuth(fixedPoint_azimuthIndex(j));        
        fixedPoint_range(i,j) = range(fixedPoint_rangeIndex(j));
        fixedPoint_dbz(i,j) = dbz_tot(fixedPoint_rangeIndex(j),fixedPoint_azimuthIndex(j));
        fixedPoint_vel(i,j) = vel(fixedPoint_rangeIndex(j),fixedPoint_azimuthIndex(j));
        fixedPoint_rhohv(i,j) = rhohv(fixedPoint_rangeIndex(j),fixedPoint_azimuthIndex(j));

    end
    
   
end

fixedPoint_Azimuth = mean(fixedPoint_azimuth);
fixedPoint_Range = mean(fixedPoint_range);


