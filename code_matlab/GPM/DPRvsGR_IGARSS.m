%DPR and GR comparison
%filename = '/home/idariash/Desktop/Denali/tmp/Ivan/Colombia/GPM_Comparison/AND/GRtoDPR.KETB.180119.022122.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/Colombia/GPM_Comparison/DROPS/GRtoDPR.KETB.180119.022122.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/PuertoRico/GPM_Comparison/GRtoDPR.TJUA.170429.017994.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/PuertoRico/GPM_Comparison/GRtoDPR.TJUA.170804.019500.V05A.DPR.NS.1_21.nc';
%filename = '/home/idariash/Desktop/Link to Ivan/PuertoRico/GPM_Comparison/GRtoDPR.TJUA.170429.017994.V05A.DPR.NS.1_21.nc';
filename = '/home/idariash/Desktop/Link to Ivan/Field_Campaigns/Relampago/Analisis/GPM/20181206/GRtoDPR.RMA1.181206.027108.V05A.DPR.NS.1_21.nc';
GR_Z = ncread(filename, 'GR_Z');
GR_Z(:,1:3) = nan;
ZFactorMeasured  = ncread(filename, 'ZFactorMeasured');
ZFactorCorrected = ncread(filename, 'ZFactorCorrected');
[r,c] = size(ZFactorMeasured);
sweep = 5;

Z_measured_DPR = reshape(ZFactorMeasured, 1, r*c);
Z_corrected_DPR = reshape(ZFactorCorrected, 1, r*c);
Z_GR = reshape(GR_Z, 1, r*c);

Z_GR(Z_GR < 18) = nan; Z_corrected_DPR(Z_GR < 18) = nan; Z_corrected_DPR(Z_GR < 18) = nan;
Z_measured_DPR(Z_measured_DPR < 18) = nan; Z_GR(Z_measured_DPR < 18) = nan;
Z_corrected_DPR(Z_corrected_DPR < 18) = nan; Z_GR(Z_corrected_DPR < 18) = nan;
% scatter(Z_GR, Z_measured_DPR)
% figure
% scatter(Z_GR, Z_corrected_DPR)

Z_GR_1st_sweep = GR_Z(:,sweep);
Z_corrected_DPR_1st_sweep = ZFactorCorrected(:,sweep);
Z_GR_1st_sweep(Z_GR_1st_sweep < 18) = nan; Z_corrected_DPR_1st_sweep(Z_GR_1st_sweep < 18) = nan;
Z_corrected_DPR_1st_sweep(Z_corrected_DPR_1st_sweep < 18) = nan; Z_GR_1st_sweep(Z_corrected_DPR_1st_sweep < 18) = nan; 

figure
scatter(Z_GR_1st_sweep, Z_corrected_DPR_1st_sweep);
x = 15:1:55;
hold on 
plot(x,x, 'r');
xlabel('Ground-based Radar Reflectivity (dBZ)');
ylabel('DPR-Ku Reflectivity (dBZ)');
grid on
CORR = nanmean((Z_GR_1st_sweep - nanmean(Z_GR_1st_sweep)).*(Z_corrected_DPR_1st_sweep - nanmean(Z_corrected_DPR_1st_sweep)))/(nanstd(Z_GR_1st_sweep)*nanstd(Z_corrected_DPR_1st_sweep));
BIAS = nanmean(Z_GR_1st_sweep - Z_corrected_DPR_1st_sweep);
RMSE = sqrt(nanmean((Z_GR_1st_sweep - Z_corrected_DPR_1st_sweep).^2));
str = ['BIAS: ' num2str(round(BIAS,2)) ' | CORR: ' num2str(round(CORR,2)) ' | RMSE: ' num2str(round(RMSE,2))];
title(str);
hold off

figure
scatter(Z_GR, Z_corrected_DPR, '.');
x = 15:1:55;
hold on 
plot(x,x, 'r');
xlabel('Ground-based Radar Reflectivity (dBZ)');
ylabel('DPR-Ku Reflectivity (dBZ)');
grid on
CORR = nanmean((Z_GR - nanmean(Z_GR)).*(Z_corrected_DPR - nanmean(Z_corrected_DPR)))/(nanstd(Z_GR)*nanstd(Z_corrected_DPR));
BIAS = nanmean(Z_GR - Z_corrected_DPR);
RMSE = sqrt(nanmean((Z_GR - Z_corrected_DPR).^2));
str = ['BIAS: ' num2str(round(BIAS,2)) ' | CORR: ' num2str(round(CORR,2)) ' | RMSE: ' num2str(round(RMSE,2))];
title(str);
hold off

%errorbar
 X = Z_GR; X(X>40) = nan;
 Y = Z_corrected_DPR; Y(Y>40) = nan;
inter_Zh=min(X):1:max(X);
%inter_Zh=10:2:34;
for i=2:length(inter_Zh)
    ss=find(X>inter_Zh(i-1) & X<inter_Zh(i) );
    x_mean(i)=(inter_Zh(i-1)+inter_Zh(i))/2;
    int_n(i)=length(ss);
    y_mean(i)=nanmean(Y(ss));
    y_std(i)=nanstd(Y(ss));
end
hold on; 
%barras = errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end),'-kx','linewidth',1.5); grid on;
barras = errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end), '-ko', 'linewidth',1.5); grid on;
%xlim([15,40])
hold off



% 
% %--------Dm------------
% GR_Dm = ncread(filename, 'GR_Dm');
% DPR_Dm = ncread(filename, 'Dm');
% GR_Dm(GR_Dm < 0) = nan;
% DPR_Dm(DPR_Dm < 0) = nan;
% x = 0:0.1:4;
% GR_Dm_1st_sweep = GR_Dm(:,sweep);
% DPR_Dm_1st_sweep = DPR_Dm(:,sweep);
% figure
% scatter(GR_Dm(:,sweep), DPR_Dm(:,sweep))
% hold on 
% plot(x,x,'r')
% xlabel('Ground-based Radar Dm (mm)');
% ylabel('DPR-Ka Dm (mm)');
% grid on
% hold off
% CORR = nanmean((GR_Dm_1st_sweep - nanmean(GR_Dm_1st_sweep)).*(DPR_Dm_1st_sweep - nanmean(DPR_Dm_1st_sweep)))/(nanstd(GR_Dm_1st_sweep)*nanstd(DPR_Dm_1st_sweep));
% BIAS = nanmean(GR_Dm_1st_sweep - DPR_Dm_1st_sweep);
% RMSE = sqrt(nanmean((GR_Dm_1st_sweep - DPR_Dm_1st_sweep).^2));
% str = ['BIAS: ' num2str(BIAS) ' | CORR: ' num2str(CORR) ' | RMSE: ' num2str(RMSE)];
% title(str);
% hold off
% Xlim = linspace(0,3,32);
% figure
% GR_Dm_hist = histogram(GR_Dm(:,sweep),Xlim);
% xlabel('Ground-based Radar Dm (mm)');
% ylabel('Counts');
% title('2017/04/29  |  Elv 1.3 deg')
% grid on
% axis([0.5 3.5 0 18000])
% 
% figure
% histogram(DPR_Dm, Xlim)
% xlabel('DPR-Ka Dm (mm)');
% ylabel('Counts');
% title('2017/04/29  |  Elv 1.3 deg')
% grid on 
% axis([0.5 3.5 0 120])
% 
% %------RR-----------
% DPR_RR = ncread(filename, 'PrecipRate');
% GR_RR_DROPS = ncread(filename, 'GR_RR_rainrate');
% DPR_RR(DPR_RR < 0) = nan;
% GR_RR_DROPS(GR_RR_DROPS < 0) = nan;
% x = 0:0.1:50;
% figure
% scatter(GR_RR_DROPS(:,sweep), DPR_RR(:,sweep))
% hold on 
% plot(x,x,'r')
% xlabel('Ground-based Radar DROPS RR (mm/h)');
% ylabel('DPR-Ka RR (mm/h)');
% grid on
% 
% GR_RR_1st_sweep = GR_RR_DROPS(:,sweep);
% DPR_RR_1st_sweep = DPR_RR(:,sweep);
% 
% CORR = nanmean((GR_RR_1st_sweep - nanmean(GR_RR_1st_sweep)).*(DPR_RR_1st_sweep - nanmean(DPR_RR_1st_sweep)))/(nanstd(GR_RR_1st_sweep)*nanstd(DPR_RR_1st_sweep));
% BIAS = nanmean(GR_RR_1st_sweep - DPR_RR_1st_sweep);
% RMSE = sqrt(nanmean((GR_RR_1st_sweep - DPR_RR_1st_sweep).^2));
% str = ['BIAS: ' num2str(BIAS) ' | CORR: ' num2str(CORR) ' | RMSE: ' num2str(RMSE)];
% title(str);
% hold off
