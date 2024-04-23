%Ivan Arias
%RELAMPAGO CHIVO
%BIRDBATH Histogram

filename = '/net/k2/storage/people/idariash/home/Field_Campaigns/PRECIP/Zdr_calibration/NetCDF/20210730/cfrad.20210730_224420.397_to_20210730_224440.397_col-radar_VER.nc';
'/net/denali/storage2/radar2/tmp/RELAMPAGO/NetCDF/20190125/cfrad.20190125_223900.602_to_20190125_223920.602_col-radar_BIRDBATH_SUR.nc';

ZDR = ncread(filename, 'ZDR');
[c , r] = size(ZDR);
ZDR = reshape(ZDR,1,r*c);
ZDR(abs(ZDR) > 2) = nan;
figure
hist(ZDR, 32)
