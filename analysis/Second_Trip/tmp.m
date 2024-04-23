% Ivan Arias
% 2020/11/09

filename_csapr = '/net/denali/storage/radar/CSAPR/DROPS/corcsapr2cfrppiqcM1.b1.20181214.013003.nc';

filename_chivo = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2018/12/14/chivo.1b.20181214_013042.REL_PFAR360.nc';

csapr_lat = ncread(filename_csapr, 'latitude');
csapr_lon = ncread(filename_csapr, 'longitude');

chivo_lat = ncread(filename_chivo, 'latitude');
chivo_lon = ncread(filename_chivo, 'longitude');

grid_lat = (chivo_lat + csapr_lat)/2;
grid_lon = (chivo_lon + csapr_lon)/2;
