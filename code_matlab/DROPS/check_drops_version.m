%2020/06/05
% Ivan Arias
% Check DROPS with different version

filename_ForDROPS_2016 = '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b_2/2019/01/13/chivo.1b.20190113_150052.PPINRLT135A.nc';
filename_ForDROPS_2020 =  '/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b_whitney/2019/01/13/chivo.1b.20190113_150052.PPINRLT135A.nc';
%'/net/denali/storage/radar/RELAMPAGO/quality_controlled_data/level_1b/2019/01/13/chivo.1b.20190113_150052.PPINRLT135A.nc';


reflectivity_2016 = ncread(filename_ForDROPS_2016, 'corrected_reflectivity');
reflectivity_2020 = ncread(filename_ForDROPS_2020, 'corrected_reflectivity');

[c, r] = size(reflectivity_2016);
reflectivity_2016 = reshape(reflectivity_2016, r*c, 1);
reflectivity_2020 = reshape(reflectivity_2020, r*c, 1);

reflectivity_2016(reflectivity_2016 < -100) = nan;

scatter(reflectivity_2016, reflectivity_2020, '.')
grid on
xlabel('Reflectivity DROPS 2016')
ylabel('Reflectivity DROPS 2020')
xlim([-40,70])
ylim([-40,70])
hold on 
x = -40:70;
plot(x,x)