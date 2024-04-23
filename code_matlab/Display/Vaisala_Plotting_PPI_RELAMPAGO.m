% Ivan Arias
% Plot Vaisala Radar products
% Nov-25/2017
% Colorado State University/ Radar and Communcation Group

addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/PPI_Plotting');
addpath('/net/denali/storage2/radar2/tmp/Ivan/Utiles/Matlab');
%file_path = '/home/idariash/Desktop/Link to Ivan/Colombia/Analisis/TAB/20180805/NetCDF/';
file_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Data/DROPS/20181126_07to08UTC/';
figure_path = '/net/denali/storage2/radar2/tmp/Ivan/Field_Campaigns/Relampago/Analisis/Self_Consistency/Figures/20181126_KDP/';
Directory = dir([file_path '/*COW*']);

rhoHV_threshold = 0.95;

for i = 1:length(Directory)
    close all
    filename = [file_path, Directory(i).name];
    RadarType = 'VAISALA';
    Campaign = 'Colombia';
    close all

    % get the fields
    Z = ncread(filename, 'DBZ');
    ZDR = ncread(filename, 'ZDR');
    VEL = ncread(filename, 'VEL');
    RHOHV = ncread(filename, 'RHOHV');
    PHIDP = ncread(filename, 'PHIDP');
    KDP = ncread(filename, 'KDP');

    % Quality filtering (SQI and rho_hv)
    Z(RHOHV < rhoHV_threshold) = NaN;
    ZDR(RHOHV < rhoHV_threshold) = NaN;
    VEL(RHOHV < rhoHV_threshold) = NaN;
    PHIDP(RHOHV < rhoHV_threshold) = NaN;
    KDP(RHOHV < rhoHV_threshold) = NaN;
    
    KDP_indexing = (KDP_rvp9 < 0.009 | KDP_rvp9 > 0.021);
    
    


    %Note: need to change the range in km not in m
    range = ncread(filename, 'range')/1e3;

    % get the elevation information
    elv = ncread(filename, 'elevation');
    azi = ncread(filename, 'azimuth');


    %configuration to show 1st sweep (lowest elevation) 360
    elv = elv(361:720);
    azi = azi(361:720);
    Z_1st_sweep = Z(360*length(range)+1:720*length(range));
    Z_1st_sweep = vec2mat(Z_1st_sweep, length(range));
    Z_1st_sweep = Z_1st_sweep';
    %Z_1st_sweep(Z_1st_sweep<0) = NaN;

    ZDR_1st_sweep = ZDR(360*length(range)+1:720*length(range));
    ZDR_1st_sweep = vec2mat(ZDR_1st_sweep, length(range));
    ZDR_1st_sweep = ZDR_1st_sweep';

    VEL_1st_sweep = VEL(360*length(range)+1:720*length(range));
    VEL_1st_sweep = vec2mat(VEL_1st_sweep, length(range));
    VEL_1st_sweep = VEL_1st_sweep';

    RHOHV_1st_sweep = RHOHV(1:360*length(range));
    RHOHV_1st_sweep = vec2mat(RHOHV_1st_sweep, length(range));
    RHOHV_1st_sweep = RHOHV_1st_sweep';

    PHIDP_1st_sweep = PHIDP(360*length(range)+1:720*length(range));
    PHIDP_1st_sweep = vec2mat(PHIDP_1st_sweep, length(range));
    PHIDP_1st_sweep = PHIDP_1st_sweep';

    KDP_1st_sweep = KDP(360*length(range)+1:720*length(range));
    KDP_1st_sweep = vec2mat(KDP_1st_sweep, length(range));
    KDP_1st_sweep = KDP_1st_sweep';


    %To find data drop via discontinuity in azi[n]
    azi(mod(diff(azi),360)>2) = nan;

    % Plot
    figure(1); %DBZ
    % aa=plotPPI(range, azi', (data.Z), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
    aa=plotPPI(range', azi', Z_1st_sweep, 'dBZ', [-30 70], .1, 50,'nwsz',1);
    [date, hour] = get_time(Directory(i).name, RadarType, Campaign);
    str = [hour, ' UTC, DBZ, ', ' rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(aa,fullfile([figure_path, '/DBZ/'], [date, '_', hour, 'dBZ_','rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(2); %ZDR 
    % aa=plotPPI(range, azi', (data.Z), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
    bb=plotPPI(range', azi', ZDR_1st_sweep, 'dBZ', [-4 4], .1, 50,'nwsz',1);
    str = [hour, ' UTC, ZDR',', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(bb,fullfile([figure_path, '/ZDR/'], [date, '_', hour, '_ZDR', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(3); %VEL
    % aa=plotPPI(range, azi', (data.Z), 'dBZ', [-30 75], 2.5, 10,'nwsz',1);
    cc=plotPPI(range', azi', VEL_1st_sweep, 'm/s', [-20 20], .01, 50,'nwsv',1);
    str = [hour, ' UTC, VEL',', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(cc,fullfile([figure_path, '/VEL/'], [date, '_', hour, 'Vel', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(4);  % RHOHV
    dd=plotPPI(range', azi', RHOHV_1st_sweep, '\rho_{hv}', [0 1], .01, 50,'jet',1);
    str = [hour, ' UTC, RhoHV'];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(dd,fullfile([figure_path, '/RHOHV/'], [date, '_', hour, 'RhoHV']),'png');

    figure(5); % PHIDP
    ee=plotPPI(range', azi', PHIDP_1st_sweep, 'deg', [-180 180], .1, 50,'jet',1);
    str = [hour, ' UTC, PHIDP', ', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    saveas(ee,fullfile([figure_path, '/PHIDP/'], [date, '_', hour, 'PHIDP', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

    figure(6); %KDP
    ff = plotPPI(range', azi', KDP_1st_sweep, 'deg/km', [-1 1], .01, 50,'winter',1);
    str = [hour, ' UTC, KDP', ', rhoHV: ', num2str(rhoHV_threshold)];
    title(str, 'Color', 'k'); 
    %grid on; 
    set(gcf, 'Position', get(0,'Screensize'));
    xlim([-150 150]);
    ylim([-150 150]);
    colormap parula
    saveas(ff,fullfile([figure_path, '/KDP/'], [date, '_', hour, 'KDP', 'rhoHV', num2str(rhoHV_threshold), '.png']),'png');

end