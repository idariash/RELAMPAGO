% Data Quality, calculate azimuth different using fixed points azimuth
% 2019/07/29 Ivan Arias

data_path = '/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20190131';
RFdifference = 14; % it should be an integer, take it from tables
directory = dir([data_path '/*HYDRO*']);


% fixedPoint_azimuthIndex = [285, 191, 239, 227, 299]; %150m range
% fixedPoint_rangeIndex = [122, 108, 99, 42, 11]; % 150m range
fixedPoint_azimuthIndex = [285, 191, 239, 227, 296]; %200m range
fixedPoint_rangeIndex = [92, 81, 74, 32, 7]; % 200m range
fixedPoint_BaseLine_azimuth = [284.0316, 190.0258, 238.0331, 226.0307, 295.1770];


fixedPoint_range = zeros(length(directory), 5);
fixedPoint_dbz = zeros(1, 5);
%fixedPoint_vel = zeros(length(directory), 5);
fixedPoint_rhohv = zeros(1, 5);

azimuth_error = nan(length(directory), 5);
azimuth_errorJustAverage = nan(length(directory), 5);
azimuth_error_total = nan(1, length(directory));
azimuth_error_totalJustAverage = nan(length(directory) , 1);
scan_time = nan(1,length(directory));
%scan_time = nan(1, length(directory));

for i =1:length(directory)
    filename = [data_path '/' directory(i).name];
    start_datetime = ncreadatt(filename, '/', 'start_datetime');
    start_datetime = start_datetime(1:length(start_datetime) - 1);
    start_datetime = datetime(start_datetime,'InputFormat','uuuu-MM-dd''T''HH:mm:ss');
    scan_time(i) = datenum(start_datetime);
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
        fixedPoint_azimuth_detection_dbz = nan(1, 7);
        for rotation = -3:3
            k = rotation + 4; % Index from 1
            azimuth_rotationIndex = fixedPoint_azimuthIndex(j) + RFdifference + rotation;
            try
                dbz_j = dbz_tot(fixedPoint_rangeIndex(j), azimuth_rotationIndex);
                rhohv_j = rhohv(fixedPoint_rangeIndex(j), azimuth_rotationIndex);
            catch 
                continue
            end
            
            if dbz_j > 45 && rhohv_j < 0.85
                disp([num2str(j) ': ' num2str(rotation) ' Ref: ' num2str(dbz_j)])
                fixedPoint_azimuth_detection(k) = azimuth(azimuth_rotationIndex);
                fixedPoint_azimuth_detection_dbz(k) = dbz_j;
            end
        end
        
        azimuth_decision = sum(~isnan(fixedPoint_azimuth_detection)); % decision to include this data in the average
        if azimuth_decision < 3
            fixedPoint_dbz = fixedPoint_azimuth_detection_dbz(~isnan(...
                fixedPoint_azimuth_detection_dbz));
            
            fixedPoint_azimuth = fixedPoint_azimuth_detection(~isnan(...
                fixedPoint_azimuth_detection));

            if azimuth_decision == 2
                fixedPoint_power = 10.^(fixedPoint_dbz/10);
                azimuth_error(i,j) = (fixedPoint_azimuth(1)*fixedPoint_power(1) + ...
                    fixedPoint_azimuth(2)*fixedPoint_power(2))/sum(fixedPoint_power) - ...
                    fixedPoint_BaseLine_azimuth(j);
            else
                azimuth_error(i,j) = nanmean(fixedPoint_azimuth_detection) - ...
                fixedPoint_BaseLine_azimuth(j);
            end
             azimuth_errorJustAverage(i,j) = nanmean(fixedPoint_azimuth_detection) - ...
                fixedPoint_BaseLine_azimuth(j);
        end
        
        
        
    end
    disp('----------')
    
    if sum(~isnan(azimuth_error(i,:)))  > 3
        azimuth_error_total(i) = nanmean(azimuth_error(i,:));
        azimuth_error_totalJustAverage(i) = ...
            nanmean(azimuth_errorJustAverage(i,:));
    end

end
azimuth_error_total = azimuth_error_total - 0.6; % Baseline bias correction
%%
dayPartition = 3;
scan_numTime = scan_time;
azimuthError_mean = nan(1, dayPartition);
azimuthError_std = nan(1, dayPartition);
azimuthError_time = nan(1, dayPartition);
for i = 1:dayPartition
    L = length(scan_time);
    azimuthError_mean(i) = nanmean(azimuth_error_total(...
        floor(L/dayPartition*(i-1)) + 1: floor(L/dayPartition*i)));
    azimuthError_std(i) = nanstd(azimuth_error_total(...
        floor(L/dayPartition*(i-1)) + 1: floor(L/dayPartition*i)));
    azimuthError_time(i) = scan_numTime(floor(L/dayPartition*(2*i-1)/2)); %median
end
figure
%scatter(scan_time, azimuth_error_total)
hold on
scatter(scan_numTime, azimuth_error_total)
errorbar(azimuthError_time, azimuthError_mean, azimuthError_std)
datetick('x', 'HH:MM')
grid on 
xlabel('HH:MM (UTC)')
ylabel('Azimuth Difference (deg.)')
title(['Drift Evolution for ' data_path(length(data_path)-7: length(data_path))])
hold off 
%%
return
%%
azimuth_correction_mean = [azimuth_correction azimuthError_mean];
azimuth_correction_std = [azimuth_correction_std azimuth azimuthError_std];
%%
errorbar(errorBar_numTime, errorBar_mean, errorBar_std)
hold on
scatter(RF_numTime, RF_error)
scatter(scan_numTime, azimuth_error_total)

errorBar_time


%azimuth_error_total = nanmean(azimuth_error, 2);

 
