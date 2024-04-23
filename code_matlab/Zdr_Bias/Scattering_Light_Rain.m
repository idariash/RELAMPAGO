%Scattering Light Rain
dataPath = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Calibration/Self_Consistency/Data/20181126/NetCDF/20181126/';
directory = dir([dataPath '/*PPICOWH360_SUR.nc']);
Z_total = nan(30, 191880);
ZDR_total = nan(30, 191880);
for i = 1:length(directory)
    filename = directory(i).name;
    Z = ncread([dataPath '/'  filename], 'DBZ');
    RHOHV =  ncread([dataPath  filename], 'RHOHV');
    range = ncread([dataPath '/' filename], 'range')/1e3;
    azimuth = ncread([dataPath  filename], 'azimuth');
    Z(RHOHV < 0.9) = nan;
    Z = Z(1:360*length(range));
    Z = vec2mat(Z, length(range))';
    Z = Z(1:533,:);
    %Z(Z < 0) = nan;
    Z(Z > 35) = nan;
    [r, c] = size(Z);
    Z = reshape(Z,1,r*c);
    Z_total(i,:) = Z;
    
    ZDR = ncread([dataPath '/' filename], 'ZDR');
    ZDR(RHOHV < 0.9) = nan;
    ZDR(ZDR<-7) = nan;
    ZDR = ZDR(1:360*length(range));
    ZDR = vec2mat(ZDR, length(range))';
    ZDR = ZDR(1:533,:);
    [r, c] = size(ZDR);
    ZDR = reshape(ZDR,1,r*c);
    %for sector computation
%     ZDR(:,azimuth<90) = nan;
%     ZDR(:,azimuth>180) = nan;
    ZDR_total(i,:) = ZDR;
end
Z = reshape(Z_total, 1, 30*191880);
ZDR = reshape(ZDR_total, 1, 30*191880);
Z=Z(~isnan(ZDR));
ZDR=ZDR(~isnan(ZDR));
