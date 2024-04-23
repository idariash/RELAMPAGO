% Test correct azimuth angle netcdf

filename = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM_CSAPR/GPM/NetCDF_tmp/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM_CSAPR/GPM/NetCDF_override/cfrad.20190113_040044.721_to_20190113_040600.628_col-radar_REL_PFAR360_SUR.nc';
azimuth = ncread(filename, 'azimuth');
azimuth_corrected = azimuth - 13.5;
ncwrite(filename, 'azimuth', azimuth_corrected);




