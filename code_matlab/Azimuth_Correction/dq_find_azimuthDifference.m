% Data Quality, calculate azimuth different using fixed points azimuth
% 2019/07/29 Ivan Arias

data_path = '/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20181129';
RFdifference = 2; % it should be an integer, take it from tables
directory = dir([data_path '/*HYDRO*']);

fixedPoint_azimuthIndex = [285, 191, 239, 227, 296];
fixedPoint_rangeIndex = [92, 81, 74, 32, 7];
fixedPoint_BaseLine_azimuth = [284.0316, 190.0258, 238.0331, 226.0307, 295.1770];


fixedPoint_range = zeros(length(directory), 5);
fixedPoint_dbz = zeros(1, 5);
%fixedPoint_vel = zeros(length(directory), 5);
fixedPoint_rhohv = zeros(1, 5);

azimuth_error = nan(length(directory), 5);
azimuth_error_total = nan(length(directory) , 1);
%scan_time = nan(1, length(directory));

for i =1:length(directory)
    filename = [data_path '/' directory(i).name];
    start_datetime = ncreadatt(filename, '/', 'start_datetime');
    start_datetime = start_datetime(1:length(start_datetime) - 1);
    scan_time(i) = datetime(start_datetime,'InputFormat','uuuu-MM-dd''T''HH:mm:ss');
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
    
    for j = 1:5
        fixedPoint_azimuth_detection = nan(1, 7);
        for rotation = -3:3
            k = rotation + 4; % Index from 1
            azimuth_rotationIndex = fixedPoint_azimuthIndex(j) + RFdifference + rotation;
            
            dbz_j = dbz_tot(fixedPoint_rangeIndex(j), azimuth_rotationIndex);
            rhohv_j = rhohv(fixedPoint_rangeIndex(j), azimuth_rotationIndex);
            
            if dbz_j > 45 && rhohv_j < 0.85
                disp([num2str(j) ': ' num2str(rotation) ' Ref: ' num2str(dbz_j)])
                fixedPoint_azimuth_detection(k) = azimuth(azimuth_rotationIndex);
                
            end
        end
        azimuth_decision = sum(~isnan(fixedPoint_azimuth_detection));
        if azimuth_decision < 3
            azimuth_error(i,j) = nanmean(fixedPoint_azimuth_detection) - ...
                fixedPoint_BaseLine_azimuth(j);
        end
    end
    disp('----------')
    
    if sum(~isnan(azimuth_error(i,:)))  > 3
        azimuth_error_total(i) = nanmean(azimuth_error(i,:));
    end

end 
%%
return
errorbar(errorBar_numTime, errorBar_mean, errorBar_std)
hold on
scatter(RF_numTime, RF_error)
scatter(scan_numTime, azimuth_error_total)


%azimuth_error_total = nanmean(azimuth_error, 2);

 
