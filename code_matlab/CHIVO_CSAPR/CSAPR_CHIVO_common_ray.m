% Ivan Arias
% 2019/09/08

% This compute the common ray between CSAPR and CHIVO

filename_chivo = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/DualRadar/Data/DROPS/cfrad.20181214_020041.798_to_20181214_020556.800_col-radar_REL_PFAR360_SUR.nc';
%filename_chivo = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/cfrad.20190125_210047.138_to_20190125_210408.475_col-radar_REL_PNL135A_PPI.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/cfrad.20190126_053009.397_to_20190126_053642.486_col-radar_PPINEARLT360_SUR_Nesbitt.nc';

filename_csapr = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/20181214/data/DROPS/corcsapr2cfrppiM1.a1.20181214.020004.nc';
%filename_csapr = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Attenuation/NetCDF/20190125_2100/DROPS/corcsapr2cfrppiM1.a1.20190125.210003.nc';
%%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/NetCDF_drops/corcsapr2cfrppiM1.a1.20190126.053003_Nesbitt.nc';

azimuth_offset = 7.4;
DBZ_chivo = ncread(filename_chivo, 'corrected_differential_reflectivity');

RHOHV_chivo = ncread(filename_chivo, 'corrected_cross_correlation_ratio');
DBZ_chivo(RHOHV_chivo < 0.8) = nan;

range_chivo = ncread(filename_chivo, 'range')'/1e3;
azimuth_chivo = ncread(filename_chivo, 'azimuth') - azimuth_offset;
elevation_chivo = ncread(filename_chivo, 'elevation');
ray_n_gates_chivo = ncread(filename_chivo, 'ray_n_gates');
delta = 0.5;
Azi_chivo = 225;
Elv_chivo = 3.3; %5; 
%6; 1.4;

Range_chivo = [];
Elevation_chivo = [];
Azimuth_chivo = [];
for I = 1:length(ray_n_gates_chivo)
    Range_chivo = [Range_chivo range_chivo(1:ray_n_gates_chivo(I))];
    
end

for I = 1:length(elevation_chivo)
    elevation_chivo_n_gates = elevation_chivo(I)*ones(1,ray_n_gates_chivo(I));
    azimuth_chivo_n_gates = azimuth_chivo(I)*ones(1,ray_n_gates_chivo(I));
    Elevation_chivo = [Elevation_chivo elevation_chivo_n_gates];
    Azimuth_chivo = [Azimuth_chivo azimuth_chivo_n_gates];
end


Index = abs(Azimuth_chivo - Azi_chivo) < delta & abs(Elevation_chivo - Elv_chivo) < delta;

dbz_chivo = DBZ_chivo(Index);
Range_chivo = Range_chivo(Index);
Elevation_chivo = Elevation_chivo(Index);
[Range_chivo, I] = sort(Range_chivo);
dbz_chivo = dbz_chivo(I);

dbz_chivo(abs(dbz_chivo) > 900) =nan;
x_projection_chivo = Range_chivo.*cosd(2.3);

% figure
% plot(x_projection_chivo, dbz_chivo)



% for CSAPR

DBZ_csapr = ncread(filename_csapr, 'corrected_differential_reflectivity');

RHOHV_csapr = ncread(filename_csapr, 'corrected_cross_correlation_ratio');
DBZ_csapr(RHOHV_csapr < 0.8) = nan;

range_csapr = ncread(filename_csapr, 'range')'/1e3;
azimuth_csapr = ncread(filename_csapr, 'azimuth');
elevation_csapr = ncread(filename_csapr, 'elevation');
ray_n_gates_csapr = ncread(filename_csapr, 'ray_n_gates');
delta = 0.5;
Azi_csapr = 45;
Elv_csapr = 2.5; %3.8;
%0.48;

Range_csapr = [];
Elevation_csapr = [];
Azimuth_csapr = [];
for I = 1:length(ray_n_gates_csapr)
    Range_csapr = [Range_csapr range_csapr(1:ray_n_gates_csapr(I))];
    
end

for I = 1:length(elevation_csapr)
    elevation_csapr_n_gates = elevation_csapr(I)*ones(1,ray_n_gates_csapr(I));
    azimuth_csapr_n_gates = azimuth_csapr(I)*ones(1,ray_n_gates_csapr(I));
    Elevation_csapr = [Elevation_csapr elevation_csapr_n_gates];
    Azimuth_csapr = [Azimuth_csapr azimuth_csapr_n_gates];
end


x_projection_chivo_commonVolume = x_projection_chivo(32 < x_projection_chivo & x_projection_chivo < 48);

Index = abs(Azimuth_csapr - Azi_csapr) < delta & abs(Elevation_csapr - Elv_csapr) < delta;

dbz_csapr = DBZ_csapr(Index);
Range_csapr = Range_csapr(Index);
Elevation_csapr = Elevation_csapr(Index);
[Range_csapr, I] = sort(Range_csapr);
dbz_csapr = dbz_csapr(I);

dbz_csapr(abs(dbz_csapr) > 900) =nan;
%%
x_projection_csapr = 77.5 - Range_csapr.*cosd(Elevation_csapr);

% % For phidp
% dbz_chivo = dbz_chivo - nanmean(dbz_chivo);
% dbz_csapr = -(dbz_csapr - nanmean(dbz_csapr));


figure
plot(x_projection_chivo, dbz_chivo)
hold on 
plot(x_projection_csapr, dbz_csapr)
hold off

grid on
xlabel('Distance from CHIVO (km)')
title('Common ray for CHIVO and CSAPR')
legend('CHIVO', 'CSAPR')
xlim([0 80])

%%
dbz_chivo_commonVolume = dbz_chivo(32 < x_projection_chivo & x_projection_chivo < 48 );
x_projection_chivo_commonVolume = x_projection_chivo(32 < x_projection_chivo & x_projection_chivo < 48);

dbz_csapr_commonVolume = dbz_csapr(32 < x_projection_csapr & x_projection_csapr < 48 );
x_projection_csapr_commonVolume = x_projection_csapr(32 < x_projection_csapr & x_projection_csapr < 48);

I = 1;
for distance_from_chivo = 32:0.2:48
    dbz_chivo_500mMean(I) = nanmean(dbz_chivo(distance_from_chivo < x_projection_chivo &...
        x_projection_chivo < distance_from_chivo + 0.2));
    dbz_csapr_500mMean(I) = nanmean(dbz_csapr(distance_from_chivo < x_projection_csapr &...
        x_projection_csapr < distance_from_chivo + 0.2));
    I = I + 1;
end

x = 25:40;

% figure
% scatter(dbz_chivo_500mMean, dbz_csapr_500mMean)
% hold on 
% plot(x,x)
% hold off
% xlim([27 37])
% ylim([27 37])

Bias = nanmean(dbz_chivo_500mMean - dbz_csapr_500mMean);
RMSE = sqrt(nanmean((dbz_chivo_500mMean - dbz_csapr_500mMean).^2));



