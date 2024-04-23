% Ivan Arias
% 2019/08/27
% Self consistency for level 1a data of CHIVO
% 

filename = '/net/denali/storage2/radar2/tmp/RELAMPAGO/quality_controlled_data/level_1a/2019/01/26/chivo.1a.20190126_053023.PPINEARLT360.nc';

%%
rhoHV_threshold = 0.98;

Z = ncread(filename, 'reflectivity');
ZDR = ncread(filename, 'differential_reflectivity');
RHOHV = ncread(filename, 'cross_correlation_ratio');
PHIDP = ncread(filename, 'differential_phase');
KDP = ncread(filename, 'specific_differential_phase');
range = ncread(filename, 'range');
elevation = ncread(filename, 'elevation');

%% get 2nd sweep

Z = Z(:, 1 < elevation & elevation < 2);
ZDR = ZDR(:,1 < elevation & elevation < 2);
RHOHV = RHOHV(:,1 < elevation & elevation < 2);
PHIDP = PHIDP(:,1 < elevation & elevation < 2);
KDP = KDP(:,1 < elevation & elevation < 2);

%% Estimation

Z_H = 10.^((Z)/10);
C = 1.46e-4;
alpha = 0.98;
beta = 0.2;
KDP_Estimated = C*(Z_H.^alpha).*10.^(-beta*ZDR)*180;
KDP_Estimated(KDP_Estimated<0.01) = nan;
KDP_Estimated(KDP_Estimated>5) = nan;
KDP(KDP<0) = nan;
Z_H_estimated =(1/C*KDP./(10.^(-beta*ZDR))).^(1/alpha);
Z_estimated = 10*log10(Z_H_estimated);

Z_estimated = reshape(Z_estimated, 993*360, 1);
Z = reshape(Z, 993*360, 1);

Z(Z<30) = nan;
Z_estimated(Z<30) = nan;
Z_estimated(Z_estimated < 30) = nan;

Z(Z>100) = nan;
Z_estimated(Z_estimated > 100) = nan;

figure
scatter(Z, Z_estimated, '.')
x = 1:70;
hold on 
plot(x,x)
xlim([30,70])
ylim([30,70])



%%
return 

close all
hist(KDP_Estimated,100)
figure
hist(KDP,100)
figure
Z_H_estimated =(1/C*KDP./(10.^(-beta*ZDR))).^(1/alpha);
Z_estimated = 10*log10(Z_H_estimated);

%%
scatter(KDP, KDP_Estimated);

