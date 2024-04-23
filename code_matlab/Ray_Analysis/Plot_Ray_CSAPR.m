% Ivan Arias
% Mayo 10/2019
% Ray over the strom 

filename = '/net/denali/storage2/radar2/tmp/Ivan/DOE/RELAMPAGO/1930_case/corcsapr2cfrhsrhiM1.a1.20190125.193715.nc';
attenuation_corrected_reflectivity_h = ncread(filename, 'attenuation_corrected_reflectivity_h');
uncorrected_reflectivity_h = ncread(filename, 'uncorrected_reflectivity_h');
range = ncread(filename, 'range')/1e3;
azimuth = ncread(filename, 'azimuth');
elevation = ncread(filename, 'elevation');
delta = 2;

attenuation_corrected_reflectivity_h_ray = attenuation_corrected_reflectivity_h(...
    :,abs(azimuth - 30) < delta & abs(elevation - 12) < delta);
attenuation_corrected_reflectivity_h_ray = nanmean(attenuation_corrected_reflectivity_h_ray,2); 

uncorrected_reflectivity_h_ray = uncorrected_reflectivity_h(...
    :,abs(azimuth - 30) < delta & abs(elevation - 12) < delta);
uncorrected_reflectivity_h_ray = nanmean(uncorrected_reflectivity_h_ray,2); 

plot(range,attenuation_corrected_reflectivity_h_ray)
hold on 
plot(range,uncorrected_reflectivity_h_ray)
hold off
