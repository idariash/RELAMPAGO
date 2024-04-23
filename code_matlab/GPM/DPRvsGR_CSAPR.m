%DPR and GR comparison
%filename = '/home/idariash/Desktop/Denali/tmp/Ivan/Colombia/GPM_Comparison/AND/GRtoDPR.KETB.180119.022122.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/Colombia/GPM_Comparison/DROPS/GRtoDPR.KETB.180119.022122.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/PuertoRico/GPM_Comparison/GRtoDPR.TJUA.170429.017994.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/PuertoRico/GPM_Comparison/GRtoDPR.TJUA.170429.017994.V05A.DPR.NS.1_21_DROPS.nc';
%close all
clear all

filename = '/net/nasstore//students/GRAD/ECE/idariash/home/GPM/2ADPR/NS/V05A/1_21/2019/GRtoDPR.CHIVO.190113.027699.V05A.DPR.NS.1_21_1b.nc';

%'/net/nasstore/students/GRAD/ECE/idariash/home/GPM/2ADPR/NS/V05A/1_21/2019/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21_b1.nc'; 

%GRtoDPR.CSAPR.190113.027699.V05A.DPR.NS.1_21.nc';
%GRtoDPR.RELAMPAGO.190113.027699.V05A.DPR.NS.1_21.nc'; 
% GRtoDPR.CSAPR.190113.027699.V05A.DPR.NS.1_21.nc';
%'/net/nasstore/students/GRAD/ECE/idariash/home/GPM/2ADPR/NS/V05A/1_21/2019/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21.nc'; 
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21_bookCoeff.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/GRtoDPR.CSAPR.181214.027231.V05A.DPR.NS.1_21.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21_bookCoeff.nc';
%'/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/GPM/CSAPR/GRtoDPR.CSAPR.190131.027991.V05A.DPR.NS.1_21.nc';
%'/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/GPM/GRtoDPR.RELAMPAGO.181206.027108.V05A.DPR.NS.1_21.nc';
sweep = 3:16;
GR_Z = ncread(filename, 'GR_Z');
GR_Z(GR_Z > 70) = nan;
GR_Z(GR_Z < 18) = nan;
%ZFactorMeasured  = ncread(filename, 'ZFactorMeasured');

ZFactorCorrected = ncread(filename, 'ZFactorCorrected');
ZFactorCorrected(ZFactorCorrected<20) = nan;

[r,c] = size(ZFactorCorrected);

%Z_measured_DPR = reshape(ZFactorMeasured, 1, r*c);
Z_corrected_DPR = reshape(ZFactorCorrected, 1, r*c);


Z_GR = reshape(GR_Z, 1, r*c);
% Z_GR(Z_GR < 20) = nan; Z_corrected_DPR(Z_GR < 20) = nan; Z_corrected_DPR(Z_GR < 20) = nan;
% %Z_measured_DPR(Z_measured_DPR < 18) = nan; Z_GR(Z_measured_DPR < 18) = nan;
% Z_corrected_DPR(Z_corrected_DPR < 20) = nan; Z_GR(Z_corrected_DPR < 20) = nan;
% scatter(Z_GR, Z_measured_DPR)
% figure
% scatter(Z_GR, Z_corrected_DPR)

Z_GR_1st_sweep = GR_Z(:,sweep);
Z_corrected_DPR_1st_sweep = ZFactorCorrected(:,sweep);
Z_GR_1st_sweep(Z_GR_1st_sweep < 18) = nan; Z_corrected_DPR_1st_sweep(Z_GR_1st_sweep < 18) = nan;
Z_corrected_DPR_1st_sweep(Z_corrected_DPR_1st_sweep < 18) = nan; Z_GR_1st_sweep(Z_corrected_DPR_1st_sweep < 18) = nan; 

[r, c] = size(Z_corrected_DPR_1st_sweep);
Z_corrected_DPR_1st_sweep = reshape(Z_corrected_DPR_1st_sweep, 1, r*c);
Z_GR_1st_sweep = reshape(Z_GR_1st_sweep, 1, r*c);


figure(2)
scatter( Z_corrected_DPR_1st_sweep, Z_GR_1st_sweep, '.');
x = 15:1:55;
hold on 
plot(x,x, 'r');
ylabel('CSU-CHIVO Reflectivity (dBZ)');
xlabel('DPR-Ku Reflectivity (dBZ)');
grid on
rho = nanmean((Z_GR_1st_sweep - nanmean(Z_GR_1st_sweep)).*(Z_corrected_DPR_1st_sweep - nanmean(Z_corrected_DPR_1st_sweep)))/(nanstd(Z_GR_1st_sweep)*nanstd(Z_corrected_DPR_1st_sweep));
Bias = nanmean(Z_corrected_DPR_1st_sweep - Z_GR_1st_sweep);
RMSE = sqrt(nanmean((Z_GR_1st_sweep - Z_corrected_DPR_1st_sweep).^2));
disp(Bias);
disp(rho)
disp(RMSE)
title(['Bias: ' num2str(round(Bias,2)) 'dB | CORR: ' num2str(round(rho,2))...
    ' | RMSE: ' num2str(round(RMSE,2)) 'dB']);

%errorbar
 Y = Z_GR_1st_sweep; Y(Y>44) = nan;
 X = Z_corrected_DPR_1st_sweep; X(X>44) = nan;
%inter_Zh=min(X):1:max(X);
inter_Zh=19:1:40;
for i=2:length(inter_Zh)
    ss=find(X>inter_Zh(i-1) & X<inter_Zh(i) );
    x_mean(i)=(inter_Zh(i-1)+inter_Zh(i))/2;
    int_n(i)=length(ss);
    y_mean(i)=nanmean(Y(ss));
    y_std(i)=nanstd(Y(ss));
end
%barras = errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end),'-kx','linewidth',1.5); grid on;
barras = errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end), '-ko', 'linewidth',1.5); grid on;
%xlim([15,40])
xlim([15,50])
ylim([15,50])
hold off