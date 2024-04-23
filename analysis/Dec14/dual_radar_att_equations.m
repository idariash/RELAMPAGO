% Ivan Arias
% Oct. 21, 2021

% Code to compute dual radar attenuation equations

%% Equations, here we go:

n = 1;
ref_threshold = 30;

% Read data
chivo_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_att_corrected_2018-12-14T02_00.nc';
chivo_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_original_2018-12-14T0200.nc';
csapr_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_att_corrected_2018-12-14T02_00.nc';
csapr_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_original_2018-12-14T02_00.nc';

chivo_corrected_ref = ncread(chivo_corrected_file, 'corrected_reflectivity');
chivo_original_ref = ncread(chivo_original_file, 'reflectivity');
csapr_corrected_ref = ncread(csapr_corrected_file, 'corrected_reflectivity');
csapr_original_ref = ncread(csapr_original_file, 'reflectivity');



% Selecting data to process

x = ncread(chivo_corrected_file, 'x')/1000;
y = ncread(chivo_corrected_file, 'y')/1000;
z = ncread(chivo_corrected_file, 'z')/1000;
[X, Y] = meshgrid(x,y);

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 

inbetween_radar = nan(301,301,20);
for i = 1:20
    inbetween_radar(:,:,i) = inbetween_radars_region;
end

chivo_corrected_ref(~inbetween_radar) = nan;
chivo_original_ref(~inbetween_radar) = nan;
csapr_corrected_ref(~inbetween_radar) = nan;
csapr_original_ref(~inbetween_radar) = nan;


chivo_corrected_ref(chivo_corrected_ref < ref_threshold) = nan;
chivo_original_ref(chivo_original_ref < ref_threshold) = nan;
csapr_corrected_ref(csapr_corrected_ref < ref_threshold) = nan;
csapr_original_ref(csapr_original_ref < ref_threshold) = nan;

chivo_attenuation = chivo_corrected_ref - chivo_original_ref; 
csapr_attenuation = csapr_corrected_ref - csapr_original_ref;
% chivo_attenuation(chivo_attenuation < -0) = nan;
% csapr_attenuation(csapr_attenuation < -0) = nan;

composite_reflectivity = (chivo_corrected_ref.*csapr_attenuation + csapr_corrected_ref.*chivo_attenuation)./...
   (csapr_attenuation + chivo_attenuation);

valid_data_number = sum(sum(sum(~isnan(composite_reflectivity))));
coordinate_dictionary = nan(valid_data_number, 3);
I = 1;
for i = 1:301
    for j = 1:301
        for k = 1:20
            if isnan(composite_reflectivity(i,j,k))
                continue
            end
            coordinate_dictionary(I, :) = [i, j, k];
            I = I + 1;
        end
    end
end

chivo_Xo = 0;
chivo_Yo = 0;
chivo_Zo = 0.450;

for i = 1:301
    for j = 1:301
        for k = 1:20
            if isnan(composite_reflectivity(i,j,k))
                continue
            end
            d = 0;
            Xo = X(i,j);
            Yo = Y(i,j);
            Zo = z(k);
            chivo_distance = sqrt((Xo - chivo_Xo)^2 + (Yo - chivo_Yo)^2 + (Zo - chivo_Zo)^2);
            V = ([Xo, Yo, Zo] - [chivo_Xo, chivo_Yo, chivo_Zo])/chivo_distance; % Normalized vector
            while d < chivo_distance
                v = d*V + [chivo_Xo, chivo_Yo, chivo_Zo]; % Line equation between CHIVO and the grid gate
                
                d = d + 0.1;
            end
        end
    end    
end
 
return

composite_reflectivity_n = composite_reflectivity(:,:,n);
composite_reflectivity_n(~inbetween_radars_region) = nan;



[r, c] = size(X);

% Case k = 1
alphas = zeros(1, 301*301*20);
for i= 1:r
    for j = 1:c
        if isnan(composite_reflectivity_n(i,j))
            continue
        end
        
        if chivo_attenuation < 0.5
            
        end
       
    end
    
end


%%
n = 1;
ref_threshold = 30;

chivo_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_att_corrected_2018-12-14T02_00.nc';
chivo_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_original_2018-12-14T0200.nc';
csapr_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_att_corrected_2018-12-14T02_00.nc';
csapr_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_original_2018-12-14T02_00.nc';

chivo_corrected_ref = ncread(chivo_corrected_file, 'corrected_reflectivity');
chivo_original_ref = ncread(chivo_original_file, 'reflectivity');
csapr_corrected_ref = ncread(csapr_corrected_file, 'corrected_reflectivity');
csapr_original_ref = ncread(csapr_original_file, 'reflectivity');

chivo_corrected_ref(chivo_corrected_ref < ref_threshold) = nan;
chivo_original_ref(chivo_original_ref < ref_threshold) = nan;
csapr_corrected_ref(csapr_corrected_ref < ref_threshold) = nan;
csapr_original_ref(csapr_original_ref < ref_threshold) = nan;

chivo_attenuation = chivo_corrected_ref - chivo_original_ref; 
csapr_attenuation = csapr_corrected_ref - csapr_original_ref;
% chivo_attenuation(chivo_attenuation < -0) = nan;
% csapr_attenuation(csapr_attenuation < -0) = nan;

composite_reflectivity = (chivo_corrected_ref.*csapr_attenuation + csapr_corrected_ref.*chivo_attenuation)./...
   (csapr_attenuation + chivo_attenuation);

x = ncread(chivo_corrected_file, 'x')/1000;
y = ncread(chivo_corrected_file, 'y')/1000;
[X, Y] = meshgrid(x,y);

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 
 
composite_reflectivity_n = composite_reflectivity(:,:,n);
composite_reflectivity_n(~inbetween_radars_region) = nan;

figure
pcolor(x, y, composite_reflectivity_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([0 70])
colormap('jet')

figure
imagesc(composite_reflectivity_n')

figure
imagesc(Y)


%%

ref_threshold = 30;
n = 1;
chivo_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_att_corrected_2018-12-14T02_00.nc';
chivo_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_original_2018-12-14T0200.nc';
csapr_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_att_corrected_2018-12-14T02_00.nc';
csapr_original_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_original_2018-12-14T02_00.nc';

chivo_corrected_ref = ncread(chivo_corrected_file, 'corrected_reflectivity');
chivo_original_ref = ncread(chivo_original_file, 'reflectivity');
csapr_corrected_ref = ncread(csapr_corrected_file, 'corrected_reflectivity');
csapr_original_ref = ncread(csapr_original_file, 'reflectivity');

chivo_corrected_ref(chivo_corrected_ref < ref_threshold) = nan;
chivo_original_ref(chivo_original_ref < ref_threshold) = nan;
csapr_corrected_ref(csapr_corrected_ref < ref_threshold) = nan;
csapr_original_ref(csapr_original_ref < ref_threshold) = nan;

chivo_attenuation = chivo_corrected_ref - chivo_original_ref; 
csapr_attenuation = csapr_corrected_ref - csapr_original_ref;

attenuation_min = min(chivo_attenuation, csapr_attenuation);

x = ncread(chivo_corrected_file, 'x')/1000;
y = ncread(chivo_corrected_file, 'y')/1000;
[X, Y] = meshgrid(x,y);

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 
    
chivo_attenuation_n = chivo_attenuation(:,:,n);
csapr_attenuation_n = csapr_attenuation(:,:,n);
attenuation_min_n = attenuation_min(:,:,n);
attenuation_min_n(~inbetween_radars_region) = nan;
attenuation_min_n(attenuation_min_n > 1) = nan;

figure
%pcolor(x, y, csapr_attenuation_n')
pcolor(x, y, attenuation_min_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([-2 10])
colormap('jet')

%%

chivo_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_att_corrected_2018-12-14T02_00.nc';
csapr_corrected_file = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_att_corrected_2018-12-14T02_00.nc';


ref_chivo = ncread(chivo_corrected_file, 'corrected_reflectivity');
ref_csapr = ncread(csapr_corrected_file, 'corrected_reflectivity');
x = ncread(chivo_corrected_file, 'x')/1000;
y = ncread(chivo_corrected_file, 'y')/1000;

[X, Y] = meshgrid(x,y);

ref_chivo(ref_chivo < -100) = nan;
ref_csapr(ref_csapr < -100) = nan;

n = 4;
ref_n_chivo = ref_chivo(:,:,n);
ref_n_csapr = ref_csapr(:,:,n);

% figure
% pcolor(x, y, ref_n')
% shading flat

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 
    

ref_n_chivo(~inbetween_radars_region) = nan;
ref_n_chivo(~inbetween_radars_region) = nan;

[c, r] = size(ref_n_chivo) ;
chivo_ref_vector = reshape(ref_n_chivo, 1, r*c);
csapr_ref_vector = reshape(ref_n_csapr, 1, r*c);

figure
scatter(chivo_ref_vector, csapr_ref_vector)

disp(nanmean(chivo_ref_vector - csapr_ref_vector))


%%
filename = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/chivo_att_corrected_2018-12-14T02_00.nc';

ref = ncread(filename, 'corrected_reflectivity');
x = ncread(filename, 'x')/1000;
y = ncread(filename, 'y')/1000;

[X, Y] = meshgrid(x,y);

ref(ref < 20) = nan;

n = 1;
ref_n = ref(:,:,n);

figure
pcolor(x, y, ref_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([0 70])
colormap('jet')

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 
    

ref_n(~inbetween_radars_region) = nan;
figure
pcolor(x, y, ref_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([0 70])
colormap('jet')
%hold on
 whole

%% CSAPR Ref

filename = '/top/students/GRAD/ECE/idariash/home/CSU/RELAMPAGO/analysis/Dec14/grided_data/20181214_0200/csapr_att_corrected_2018-12-14T02_00.nc';

ref = ncread(filename, 'corrected_reflectivity');
x = ncread(filename, 'x')/1000;
y = ncread(filename, 'y')/1000;

[X, Y] = meshgrid(x,y);

ref(ref < 20) = nan;

n = 1;
ref_n = ref(:,:,n);

figure
pcolor(x, y, ref_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([0 70])
colormap('jet')

Y_csapr_upper = 2.1445*X + 61.8030;
Y_csapr_lower = 0.4663*X - 28.8198;

Y_chivo_upper = 0.4663*X;
Y_chivo_lower = 2.1450*X;

inbetween_radars_region = Y_csapr_upper > Y & Y_csapr_lower < Y & ...
    Y_chivo_upper > Y & Y_chivo_lower < Y; 
    

ref_n(~inbetween_radars_region) = nan;
figure
pcolor(x, y, ref_n')
shading flat
xlim([-60 0])
ylim([-60 0])
hc=colorbar;
caxis([0 70])
colormap('jet')
%hold on


