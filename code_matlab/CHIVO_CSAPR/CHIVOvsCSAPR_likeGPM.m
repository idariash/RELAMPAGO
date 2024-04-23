% Ivan Arias
% 2019/09/07
% Scattergram CHIVO CSAPR

addpath /net/denali/storage2/radar2/people/idariash/home/Utiles/Matlab
clear all

% Nov 30 stratiform case
%   Paper
filename_chivoGrid = '/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/grid_paper/chivo_grid_20181130-0330_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/RMA_CSAPR/grid/rma_grid_20181130-0315_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/grid_paper/chivo_grid_20181130-0330_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/RMA1/20181130/grid/chivo_grid_20181130-0400_disdrometer.nc';
filename_csaprGrid = '/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/grid_paper/csapr_grid_20181130-0330_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/RMA_CSAPR/grid/csapr_grid_20181130-0315_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/DROPS/grid_paper/csapr_grid_20181130-0330_disdrometer.nc';
%'/net/denali/storage2/radar2/people/idariash/home/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/RMA1/20181130/grid/rma_grid_20181130-0400_disdrometer.nc';

% filename_chivoGrid = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/chivo_grid_20181130-0330.nc';
% filename_csaprGrid = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/csapr_grid_20181130-0330.nc';
% to compute the grid, I use pyart code in: 

% Jan 26 stratiform case presented in AMS Radar
% filename_chivoGrid = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/chivo_grid_20190126-0530_Book.nc';
% filename_csaprGrid = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/CHIVO_CSAPR/data/csapr_grid_20190126-0530_Book.nc';

chivo_reflectivity = ncread(filename_chivoGrid, 'corrected_reflectivity');
chivo_diffRef = ncread(filename_chivoGrid, 'corrected_differential_reflectivity');
chivo_rhohv = ncread(filename_chivoGrid, 'corrected_copol_correlation_ratio');

csapr_reflectivity = ncread(filename_csaprGrid, 'corrected_reflectivity');
csapr_diffRef = ncread(filename_csaprGrid, 'corrected_differential_reflectivity');
csapr_rhohv = ncread(filename_csaprGrid, 'corrected_copol_correlation_ratio');

chivo_reflectivity(chivo_rhohv < 0.9) = nan;
csapr_reflectivity(csapr_rhohv < 0.9) = nan;
chivo_diffRef(chivo_rhohv < 0.9) = nan;
csapr_diffRef(csapr_rhohv < 0.9) = nan;

% chivo_reflectivity(chivo_reflectivity < 25) = nan;
% csapr_reflectivity(csapr_reflectivity < 25) = nan;
% chivo_diffRef(chivo_reflectivity < 25) = nan;
% csapr_diffRef(csapr_reflectivity < 25) = nan;

chivo_reflectivity = reshape(chivo_reflectivity, 51*51, 7);
chivo_DiffRef = reshape(chivo_diffRef, 51*51, 7);
chivo_rhohv = reshape(chivo_rhohv, 51*51, 7);

csapr_reflectivity = reshape(csapr_reflectivity, 51*51, 7);
csapr_DiffRef = reshape(csapr_diffRef, 51*51, 7);
csapr_rhohv = reshape(csapr_rhohv, 51*51, 7);




x = -50:50;
level = 4;

chivo_reflectivity_vector = [chivo_reflectivity(:,2); chivo_reflectivity(:,3); chivo_reflectivity(:,4)];
csapr_reflectivity_vector = [csapr_reflectivity(:,2); csapr_reflectivity(:,3); csapr_reflectivity(:,4)];

chivo_DiffRef_vector = [chivo_DiffRef(:,2); chivo_DiffRef(:,3); chivo_DiffRef(:,4)];
csapr_DiffRef_vector = [csapr_DiffRef(:,2); csapr_DiffRef(:,3); csapr_DiffRef(:,4)];

rho = nanmean((chivo_reflectivity_vector - nanmean(chivo_reflectivity_vector)).*(csapr_reflectivity_vector - nanmean(csapr_reflectivity_vector)))/(nanstd(chivo_reflectivity_vector)*nanstd(csapr_reflectivity_vector));
Bias = nanmean(chivo_reflectivity_vector - csapr_reflectivity_vector);
RMSE = sqrt(nanmean((chivo_reflectivity_vector - csapr_reflectivity_vector).^2));

figure

%errorbar
 Y = csapr_reflectivity_vector; %Y(Y>44) = nan;
 X = chivo_reflectivity_vector; %X(X>44) = nan;
%inter_Zh=min(X):1:max(X);
inter_Zh=17:4:42;
for i=2:length(inter_Zh)
    ss=find(X>inter_Zh(i-1) & X<inter_Zh(i) );
    x_mean(i)=(inter_Zh(i-1)+inter_Zh(i))/2;
    int_n(i)=length(ss);
    y_mean(i)=nanmean(Y(ss));
    y_std(i)=nanstd(Y(ss));
end
%barras = errorbar(x_mean(2:end),y_mean(2:end),y_std(2:end),'-kx','linewidth',1.5); grid on;
barras = errorbar(x_mean(3:end),y_mean(3:end),y_std(3:end), '-ko', 'linewidth',1.5); grid on;
%xlim([15,40])
xlim([15,50])
ylim([15,50])
hold on


% [A,CZplot,CZDRplot]=p_icp(chivo_reflectivity_vector, csapr_reflectivity_vector,50);
% mean_A = nanmean(A,2);
% pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
% caxis([0 2]);
% hold on 
plot(x,x, 'r')
% hc=colorbar;
% xlabel(hc,'log_{10}(N_{obs})');
xlim([20 50])
ylim([20 50])
title(['Bias: ' num2str(round(Bias,2)) 'dB | CORR: ' num2str(round(rho,2))...
    ' | RMSE: ' num2str(round(RMSE,2)) 'dB']);
xlabel('CHIVO Reflectivity (dBZ)')
ylabel('CSAPR Reflectivity (dBZ)')

return

rho_Zdr = nanmean((chivo_DiffRef_vector - nanmean(chivo_DiffRef_vector)).*(csapr_DiffRef_vector - nanmean(csapr_DiffRef_vector)))/(nanstd(chivo_DiffRef_vector)*nanstd(csapr_DiffRef_vector));
Bias_Zdr = nanmean(csapr_DiffRef_vector - chivo_DiffRef_vector);
RMSE_Zdr = sqrt(nanmean((chivo_DiffRef_vector - csapr_DiffRef_vector).^2));

[A,CZplot,CZDRplot]=p_icp(chivo_DiffRef_vector, csapr_DiffRef_vector,50);
mean_A = nanmean(A,2);
figure; pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
caxis([0 2]);
hold on 
plot(x,x, 'r')
hc=colorbar;
xlabel(hc,'log_{10}(N_{obs})');
xlim([0 2.5])
ylim([0 2.5])
title(['Bias: ' num2str(round(Bias_Zdr,2)) 'dB | CORR: ' num2str(round(rho_Zdr,2))...
    ' | RMSE: ' num2str(round(RMSE_Zdr,2)) 'dB']);
xlabel('CHIVO Diff. Reflectivity (dB)')
ylabel('CSAPR Diff. Reflectivity (dB)')

return

figure
scatter(chivo_reflectivity_vector, csapr_reflectivity_vector, '.');
hold on 
plot(x,x)
xlim([20 50])

xlabel('CHIVO Reflectivity (dBZ)');
ylabel('CSAPR Reflectivity (dBZ)');
grid on
disp(Bias);
disp(rho)
disp(RMSE)
title(['Bias: ' num2str(round(Bias,2)) 'dB | CORR: ' num2str(round(rho,2))...
    ' | RMSE: ' num2str(round(RMSE,2)) 'dB']);
xlabel('CHIVO Reflectivity (dBZ)')
ylabel('CSAPR Reflectivity (dBZ)')

%%
[A,CZplot,CZDRplot]=p_icp(chivo_reflectivity_vector, csapr_reflectivity_vector,50);
mean_A = nanmean(A,2);
figure; pcolor(CZplot, CZDRplot, log10(A)); shading flat; grid on; hold on
caxis([0 2]);
xlabel('RMA Reflectivity (dBZ)')
ylabel('CSAPR Reflectivity (dBZ)')
return 
my_handle=colorbar('ytick',1:0.2:5.4,'FontWeight','bold');
set(get(my_handle,'Title'),'string','log_{10}(N_{obs})','FontWeight','bold');
axis([10 40 -4 4]);
set(gcf, 'Position', [100, 100, 800, 800])
set(gca,'XTick',[10:2:50],'FontSize',20,'FontWeight','bold');
set(gca,'YTick',[-4:0.5:4],'FontSize',20,'FontWeight','bold');
xlabel(' Z (dBZ)','FontSize',20,'FontWeight','bold');
ylabel(' ZDR (dB)','FontSize',20,'FontWeight','bold');

% figure
% scatter(chivo_DiffRef_vector, csapr_DiffRef_vector, '.')
% hold on 
% plot(x,x)
% xlim([0 2])
% hold off
% 
% figure
% scatter(chivo_rhohv(:,level), csapr_rhohv(:,level), '.')
% hold on 
% plot(x,x)
% xlim([0.95 1])
% hold off