% Ivan Arias
% Plot Vaisala Radar products
% Nov-25/2017
% Colorado State University/ Radar and Communcation Group

addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting');
addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab');
%file_path = '/home/idariash/Desktop/Link to Ivan/Colombia/Analisis/TAB/20180805/NetCDF/';
file_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Doppler/DROPS/';
figure_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Figures/20181126_KDP/';
Directory = dir([file_path '/*RHI*']);

rhoHV_threshold = 0.8;



for i = 1:1 %length(Directory)
    close all
    filename = [file_path, Directory(i).name];
    RadarType = 'VAISALA';
    Campaign = 'Colombia';
    close all

    % get the fields
    Z = ncread(filename, 'DBZ');
    ZDR = ncread(filename, 'ZDR');
    %VEL = ncread(filename, 'VEL');
    RHOHV = ncread(filename, 'RHOHV');
    PHIDP = ncread(filename, 'PHIDP');
    KDP = ncread(filename, 'KDP');
    HyC = ncread(filename, 'HydroClass');

    % Quality filtering (SQI and rho_hv)
    Z(RHOHV < rhoHV_threshold) = NaN;
    ZDR(RHOHV < rhoHV_threshold) = NaN;
    %VEL(RHOHV < rhoHV_threshold) = NaN;
    PHIDP(RHOHV < rhoHV_threshold) = NaN;
    KDP(RHOHV < rhoHV_threshold) = NaN;
    
   % KDP_indexing = (KDP_rvp9 < 0.009 | KDP_rvp9 > 0.021);
    
    


    %Note: need to change the range in km not in m
    range = ncread(filename, 'range')/1e3;

    % get the elevation information
    elv = ncread(filename, 'elevation');
    azi = ncread(filename, 'azimuth');


    %configuration to show 1st sweep (lowest elevation) 360
    elv = elv(2309:2484); % change this depending on the azimuth
    azi = azi(2309:2484);
    Z_1st_sweep = Z(2308*length(range) + 1: 2484*length(range));
    Z_1st_sweep = vec2mat(Z_1st_sweep, length(range));
    Z_1st_sweep = Z_1st_sweep';
    %Z_1st_sweep(Z_1st_sweep<0) = NaN;

    ZDR_1st_sweep = ZDR(2308*length(range) + 1: 2484*length(range));
    ZDR_1st_sweep = vec2mat(ZDR_1st_sweep, length(range));
    ZDR_1st_sweep = ZDR_1st_sweep';

    RHOHV_1st_sweep = RHOHV(2308*length(range) + 1: 2484*length(range));
    RHOHV_1st_sweep = vec2mat(RHOHV_1st_sweep, length(range));
    RHOHV_1st_sweep = RHOHV_1st_sweep';

    PHIDP_1st_sweep = PHIDP(2308*length(range) + 1: 2484*length(range));
    PHIDP_1st_sweep = vec2mat(PHIDP_1st_sweep, length(range));
    PHIDP_1st_sweep = PHIDP_1st_sweep';

    KDP_1st_sweep = KDP(2308*length(range) + 1: 2484*length(range));
    KDP_1st_sweep = vec2mat(KDP_1st_sweep, length(range));
    KDP_1st_sweep = KDP_1st_sweep';
    
    HyC_1st_sweep = HyC(2308*length(range) + 1: 2484*length(range));
    HyC_1st_sweep = vec2mat(HyC_1st_sweep, length(range));
    HyC_1st_sweep = HyC_1st_sweep';


    %To find data drop via discontinuity in azi[n]
    %azi(mod(diff(azi),360)>2) = nan;

    % Plot
    figure(1); %DBZ
    %plotRHI2(range', elv', Z_1st_sweep, 'dBZ', [-10 80], 0.1, 50,'',1);
    plotRHI2(range', elv', Z_1st_sweep, 'dBZ', [-40 80], .1, 20, 2, 'nwsz',1);
    %plotRHI2(range', elv', double(HyC_1st_sweep), '', [-0.5 18], 1, 5,2 ,'',1);
    [date, hour] = get_time(Directory(i).name, RadarType, Campaign);
    str = [hour, ' UTC, DBZ, ', ' rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    %saveas(aa,fullfile([figure_path, '/DBZ/'], [date, '_', hour, 'dBZ_','rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(2); %ZDR 
    % aaplotRHI2(range, elv', (data.Z), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
    plotRHI2(range', elv', ZDR_1st_sweep, 'dB', [-4 8], .1, 20, 2,'jet',1);
    str = [hour, ' UTC, ZDR',', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    %saveas(bb,fullfile([figure_path, '/ZDR/'], [date, '_', hour, '_ZDR', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

%     figure(3); %VEL
%     % aaplotRHI2(range, elv', (data.Z), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
%     ccplotRHI2(range', elv', VEL_1st_sweep, 'm/s', [-20 20], .01, 50,'nwsv',1);
%     str = [hour, ' UTC, VEL',', rhoHV: ', num2str(rhoHV_threshold)];
%     title(str, 'Color', 'k'); 
%     %grid on; 
%     set(gcf, 'Position', get(0,'Screensize'));
%     xlim([0 120]);
%     ylim([0 21]);
%     saveas(cc,fullfile([figure_path, '/VEL/'], [date, '_', hour, 'Vel', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(4);  % RHOHV
    plotRHI2(range', elv', RHOHV_1st_sweep, '\rho_{hv}', [0.9 1], .1, 20, 2,'jet',1);
    str = [hour, ' UTC, RhoHV'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    %saveas(dd,fullfile([figure_path, '/RHOHV/'], [date, '_', hour, 'RhoHV']),'png');

    figure(5); % PHIDP
    plotRHI2(range', elv', PHIDP_1st_sweep, 'deg', [-180 180], .1, 20, 2,'jet',1);
    str = [hour, ' UTC, PHIDP' ];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    %saveas(ee,fullfile([figure_path, '/PHIDP/'], [date, '_', hour, 'PHIDP', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(6); %KDP
    plotRHI2(range', elv', KDP_1st_sweep, 'deg/km', [-1 8], .1, 20, 2,'winter',1);
    str = [hour, ' UTC, KDP' ];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    colormap jet
    %saveas(ff,fullfile([figure_path, '/KDP/'], [date, '_', hour, 'KDP', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');
    
    figure(7); %HyC
    HyC_1st_sweep(HyC_1st_sweep < 5) = nan;
    plotRHI2(range', elv', double(HyC_1st_sweep), '', [-0.5 18], 1, 20, 2,'',1);
    units = 'HyC';
    color_map = 'HSV';
    CAxsLim = [5 17];
    incrC = 16/16;
    hcb=fcolorbar(CAxsLim(1):incrC:CAxsLim(2),units,color_map);
    %colormap(hc_map)
    set(hcb,'YTick',floor(1:16),'FontSize',9); %min and max val
    set(hcb,'TickLength',[0 0]);
    
    colorbar('Ticks',[5.5, 6.5, 7.5, 8.5, 9.5, 10.5 11.5 12.5 13.5 14.5 15.5 16.5],...
             'TickLabels',{'ND','LD','DR','RA','HR', 'RH', 'HA', 'GR', 'WS', 'DI', 'CR', 'DN'});
    
    str = [hour, ' UTC, HydroClass' ];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([0 120]);
    ylim([0 21]);
    colormap HSV
end